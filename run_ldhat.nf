#! /usr/bin/env nextflow
nextflow.enable.dsl = 2


process LDHAT_REFORMAT_FASTA{
    // publishDir "Output", mode: "copy", saveAs: {filename -> "${fn_modifier}_${filename}"}

    input:
        path fasta

    output:
        tuple stdout,
            path("LDhat_reformated.fa")

    script:
    """
    LDhat_reformat_fasta.py ${fasta}
    """
}


process LDHAT_CONVERT{
//     publishDir "Output", mode: "copy", saveAs: {filename -> "${fn_modifier}_${filename}"}

    input:
        tuple val(sample_size),
            path("LDhat_reformated.fa")

    output:
        tuple val(sample_size),
            path("locs.txt"),
            path("sites.txt")

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
    // publishDir "Output", mode: "copy", saveAs: {filename -> "${fn_modifier}_${filename}"}

    input:
        tuple val(sample_size),
            path("locs.txt"),
            path("sites.txt")

    output:
        tuple val(sample_size),
            path("sites.txt"),
            path("locs_C.txt")

    script:
        """
        get_ldhat_to_use_gene_conv.py locs.txt
        """

}

process WATTERSON_ESTIMATE {
    // publishDir "Output", mode: "copy", saveAs: {filename -> "${fn_modifier}_${filename}"}

    input:
        tuple val(sample_size),
            path("sites.txt"),
            path("locs_C.txt")

    output:
        tuple stdout,
            val(sample_size),
            path("sites.txt"),
            path("locs_C.txt")

    script:
    """
    snps=\$(head -1 locs_C.txt | awk '{split(\$0,a," "); print a[1]}')
    genome_len=\$(head -1 locs_C.txt | awk '{split(\$0,a," "); print a[2]}' | awk '{split(\$0,a,"."); print a[1]}')
    watterson_estimate.py \$snps \$genome_len ${sample_size}
    """
}


process LOOKUP_TABLE_LDPOP {
    publishDir "Output", mode: "copy", saveAs: {filename -> "${fn_modifier}_${filename}"}

    input:
        tuple val(theta_est),
            val(sample_size),
            path("sites.txt"),
            path("locs_C.txt")

    output:
        tuple val(theta_est),
            val(sample_size),
            path("sites.txt"),
            path("locs_C.txt"),
            path("lookupTable.txt")

    script:
    // There are other parameters that can be adjusted, I've left them out for the time being
    // also they mention twice muation and recom rate, for the mutation and recom parameters which I am unsure how to interpret
    """
    ldtable.py --cores $task.cpus -n ${sample_size} -th ${theta_est} -rh ${params.ldpop_rho_range} --approx > lookupTable.txt
    """
}


process DOWNSAMPLED_LOOKUP_TABLE {
    // publishDir "Output", mode: "copy", saveAs: {filename -> "${fn_modifier}_${filename}"}

    input:
        tuple val(theta_est),
            val(sample_size),
            path("sites.txt"),
            path("locs_C.txt")
        
        val mutation_rate
        path downsampled_lookup_tables
        

    output:
        tuple val(theta_est),
            val(sample_size),
            path("sites.txt"),
            path("locs_C.txt"),
            path("lookupTable.txt")

    script:
    """
    reformat_downsampled_lk_table.py lk_downsampled_${sample_size}.csv ${sample_size} ${mutation_rate} ${params.ldpop_rho_range}
    """
}


process LDHAT_PAIRWISE{
    publishDir "Output", mode: "copy", saveAs: {filename -> "${fn_modifier}_${filename}"}

    input:
        tuple val(theta_est),
            val(sample_size),
            path("sites.txt"),
            path("locs_C.txt"),
            path("lookupTable.txt")

    output:
        tuple path("pairwise_freqs.txt"),
            path("pairwise_outfile.txt"),
            path("pairwise_stdOut.txt"),
            path("watterson_est.txt")

    script:
        // uses pexpect to handle unavoidable prompts
        """
        run_pairwise_with_pexpect.py ${params.recom_tract_len} sites.txt locs_C.txt lookupTable.txt > pairwise_stdOut.txt
        echo ${theta_est} > watterson_est.txt
        """

}


process PAIRWISE_PROCESS_OUTPUT{
    publishDir "Output", mode: "copy", saveAs: {filename -> "${fn_modifier}_${filename}"}

    input:
        tuple val(rho), 
            val(theta),
            val(genome_size),
            val(sample_size),
            val(seed), 
            val(fn_modifier),
            path("pairwise_freqs.txt"),
            path("pairwise_outfile.txt"),
            path("pairwise_stdOut.txt")

    output:
        path "processed_results.csv", emit: processed_results_csv

    script:
        """
        pairwise_process_output.py pairwise_outfile.txt ${rho} ${sample_size} ${genome_size}
        """

}

process PLOT_RESULTS{
    publishDir "Output/Results", mode: "copy"

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

    params.mutation_rate = 0.01
    params.recom_tract_len = 2300
    params.ldpop_rho_range = "101,100"

    params.input_fasta = 'none'
    params.lookup_tables = "Lookup_tables"

    // Input verification
    if (params.input_fasta == 'none') {
        println "No input .fa specified. Use --input_fasta [.fa]"
        exit 1
    }

    // Channels
    input_fasta_channel = Channel.fromPath( params.input_fasta )
    downsampled_lookup_tables = Channel.fromPath( "${params.lookup_tables}/lk_downsampled_*.csv" ).collect()

    // For each process there is a output of tuple with the params that change + necessary files/values  to move forward until they are no longer need

    LDHAT_REFORMAT_FASTA(input_fasta_channel)

    LDHAT_CONVERT(LDHAT_REFORMAT_FASTA.out)

    SWITCH_TO_GENE_CONVERSION_MODE(LDHAT_CONVERT.out)

    WATTERSON_ESTIMATE(SWITCH_TO_GENE_CONVERSION_MODE.out)

    // LOOKUP_TABLE_LDPOP(WATTERSON_ESTIMATE.out)

    DOWNSAMPLED_LOOKUP_TABLE(WATTERSON_ESTIMATE.out, params.mutation_rate, downsampled_lookup_tables)

    LDHAT_PAIRWISE(DOWNSAMPLED_LOOKUP_TABLE.out)

    // PAIRWISE_PROCESS_OUTPUT(LDHAT_PAIRWISE.out)

    // collectedFile = PAIRWISE_PROCESS_OUTPUT.out.processed_results_csv.collectFile(name:"collected_results.csv",storeDir:"Output/Results", keepHeader:true)

    // PLOT_RESULTS(collectedFile)
}