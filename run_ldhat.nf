#! /usr/bin/env nextflow
nextflow.enable.dsl = 2


process PREFIX_FILENAME {

    // maxForks 1

    // echo true

    input:
        path(fasta)
        val(prefix_fn)

    output:
        tuple stdout,
            path(fasta)

    script:
    """
    prefix_filename.py ${fasta} ${prefix_fn} 
    """
}


process LDHAT_REFORMAT_FASTA{
    // publishDir "LDhat_Output", mode: "copy", saveAs: {filename -> "${prefix_filename}${filename}"}

    input:
        tuple val(prefix_filename),
            path(fasta)

    output:
        tuple val(prefix_filename),
            stdout,
            path("LDhat_reformated.fa")

    script:
    """
    LDhat_reformat_fasta.py ${fasta}
    """
}


process LDHAT_CONVERT{
//     publishDir "LDhat_Output", mode: "copy", saveAs: {filename -> "${prefix_filename}${filename}"}

    input:
        tuple val(prefix_filename),
            val(sample_size),
            path("LDhat_reformated.fa")

    output:
        tuple val(prefix_filename),
            val(sample_size),
            path("locs.txt"),
            path("sites.txt"),
            path("LDhat_reformated.fa")

    script:
        // -2only: Specifies that only polymorphic sites with exactly two alleles
        // will be analysed and outputted Although only those sites with two alleles
        // are analysed in pairwise and interval, outputting all segregating sites
        // may be of interest and can be used to estimate a finite-sites estimate of
        // Watterson’s theta per site within pairwise
        """
        convert -seq LDhat_reformated.fa -2only
        """

}


process SWITCH_TO_GENE_CONVERSION_MODE {
    // publishDir "LDhat_Output", mode: "copy", saveAs: {filename -> "${prefix_filename}${filename}"}

    input:
        tuple val(prefix_filename),
            val(sample_size),
            path("locs.txt"),
            path("sites.txt"),
            path("LDhat_reformated.fa")

    output:
        tuple val(prefix_filename),
            val(sample_size),
            path("sites.txt"),
            path("locs_C.txt"),
            path("LDhat_reformated.fa")

    script:
        """
        get_ldhat_to_use_gene_conv.py locs.txt
        """

}

process WATTERSON_ESTIMATE {
    publishDir "LDhat_Output", mode: "copy", saveAs: {filename -> "${prefix_filename}${filename}"}

    input:
        tuple val(prefix_filename),
            val(sample_size),
            path("sites.txt"),
            path("locs_C.txt"),
            path("LDhat_reformated.fa")

    output:
        tuple val(prefix_filename),
            stdout,
            val(sample_size),
            path("sites.txt"),
            path("locs_C.txt")

    script:
    """
    fasta_variant_sites.py LDhat_reformated.fa
    snps=\$(wc -l variants_in_fasta.csv | awk '{split(\$0,a," "); print a[1]}')
    genome_len=\$(head -1 locs_C.txt | awk '{split(\$0,a," "); print a[2]}' | awk '{split(\$0,a,"."); print a[1]}')
    original_watterson.py \$genome_len \$snps ${sample_size}
    """
}


process LOOKUP_TABLE_LDPOP {
    publishDir "LDhat_Output", mode: "copy", saveAs: {filename -> "${prefix_filename}${filename}"}

    input:
        tuple val(prefix_filename),
            val(theta_est),
            val(sample_size),
            path("sites.txt"),
            path("locs_C.txt")

    output:
        tuple val(prefix_filename),
            val(theta_est),
            val(sample_size),
            path("sites.txt"),
            path("locs_C.txt"),
            path("lookupTable.txt")

    script:
    // There are other parameters that can be adjusted, I've left them out for the time being
    // also they mention twice muation and recom rate, for the mutation and recom parameters which I am unsure how to interpret
    """
    ldtable.py --cores $task.cpus -n ${sample_size} -th ${theta_est} -rh ${params.lookup_grid} --approx > lookupTable.txt
    """
}


process DOWNSAMPLED_LOOKUP_TABLE {
    // publishDir "LDhat_Output", mode: "copy", saveAs: {filename -> "${prefix_filename}${filename}"}

    input:
        tuple val(prefix_filename),
            val(theta_est),
            val(sample_size),
            path("sites.txt"),
            path("locs_C.txt")
        
        val theta
        path downsampled_lookup_tables
        

    output:
        tuple val(prefix_filename),
            val(theta_est),
            val(sample_size),
            path("sites.txt"),
            path("locs_C.txt"),
            path("lookupTable.txt")

    script:
    """
    reformat_downsampled_lk_table.py lk_downsampled_${sample_size}.csv ${sample_size} ${theta} ${params.lookup_grid}
    """
}


process LDHAT_PAIRWISE{
    publishDir "LDhat_Output", mode: "copy", saveAs: {filename -> "${prefix_filename}${filename}"}

    input:
        tuple val(prefix_filename),
            val(theta_est),
            val(sample_size),
            path("sites.txt"),
            path("locs_C.txt"),
            path("lookupTable.txt")

    output:
        tuple val(prefix_filename),
            path("pairwise_freqs.txt"),
            path("pairwise_outfile.txt"),
            path("pairwise_stdOut.txt"),
            path("theta_est.csv")

    script:
        // uses pexpect to handle unavoidable prompts
        """
        run_pairwise_with_pexpect.py ${params.recom_tract_len} sites.txt locs_C.txt lookupTable.txt > pairwise_stdOut.txt
        echo ${theta_est} > theta_est.csv
        """

}


process PAIRWISE_PROCESS_OUTPUT{
    publishDir "LDhat_Output", mode: "copy", saveAs: {filename -> "${prefix_filename}${filename}"}

    input:
        tuple val(prefix_filename),
            path("pairwise_freqs.txt"),
            path("pairwise_outfile.txt"),
            path("pairwise_stdOut.txt"),
            path("theta_est.csv")

    output:
        path "processed_results.csv", emit: processed_results_csv

    script:
        """
        pairwise_process_output.py pairwise_outfile.txt
        """

}

process PLOT_RESULTS{
    publishDir "LDhat_Output/Results", mode: "copy"

    input:
        path collectedFile

    output:
        path "rho_comparison.png"
        path "max_lk_comparison.png"

    script:
        """
        plot_results.py collected_results.csv
        """

}


workflow {
    // Note: Channels can be called unlimited number of times in DSL2
    // A process component can be invoked only once in the same workflow context

    // params.theta = 0.01 // scaled theta
    params.recom_tract_len = 1000
    params.lookup_grid = "101,100" // The range of rho values used to generate lookup tables

    params.prefix_filename = 'none' // prefix string to output filenames to help distinguish runs
    params.input_fasta = 'none'
    // params.lookup_tables = "Lookup_tables" // Location of downsampled lookup tables

    // Input verification
    if (params.input_fasta == 'none') {
        println "No input .fa specified. Use --input_fasta [.fa]"
        exit 1
    }

    // Channels
    input_fasta_channel = Channel.fromPath( params.input_fasta )
    // downsampled_lookup_tables = Channel.fromPath( "${params.lookup_tables}/lk_downsampled_*.csv" ).collect()

    // For each process there is a output of tuple with the params that change + necessary files/values  to move forward until they are no longer need
    PREFIX_FILENAME(input_fasta_channel, params.prefix_filename)

    LDHAT_REFORMAT_FASTA(PREFIX_FILENAME.out)

    LDHAT_CONVERT(LDHAT_REFORMAT_FASTA.out)

    SWITCH_TO_GENE_CONVERSION_MODE(LDHAT_CONVERT.out)

    WATTERSON_ESTIMATE(SWITCH_TO_GENE_CONVERSION_MODE.out)

    LOOKUP_TABLE_LDPOP(WATTERSON_ESTIMATE.out)

    // DOWNSAMPLED_LOOKUP_TABLE(WATTERSON_ESTIMATE.out, params.theta, downsampled_lookup_tables)

    LDHAT_PAIRWISE(LOOKUP_TABLE_LDPOP.out)

    PAIRWISE_PROCESS_OUTPUT(LDHAT_PAIRWISE.out)

    // collectedFile = PAIRWISE_PROCESS_OUTPUT.out.processed_results_csv.collectFile(name:"collected_results.csv",storeDir:"Output/Results", keepHeader:true)

    // PLOT_RESULTS(collectedFile)
}