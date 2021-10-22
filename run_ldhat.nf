#! /usr/bin/env nextflow
nextflow.enable.dsl = 2


process LDHAT_REFORMAT_FASTA{
    publishDir "Output", mode: "copy", saveAs: {filename -> "${fn_modifier}_${filename}"}

    input:
        path fasta
        val mutation_rate

    output:
        tuple val(mutation_rate),
            stdout,
            path("LDhat_reformated.fa")

    script:
    """
    LDhat_reformat_fasta.py ${fasta}
    """
}


process LOOKUP_TABLE_LDPOP {
//     publishDir "Output", mode: "copy", saveAs: {filename -> "${fn_modifier}_${filename}"}

    input:
        tuple val(mutation_rate),
            val(sample_size),
            path("LDhat_reformated.fa")

    output:
        tuple val(mutation_rate),
            val(sample_size),
            path("LDhat_reformated.fa"),
            path("lookupTable.txt")

    script:
    // There are other parameters that can be adjusted, I've left them out for the time being
    // also they mention twice muation and recom rate, for the mutation and recom parameters which I am unsure how to interpret
    """
    ldtable.py --cores 4 -n ${sample_size} -th ${mutation_rate} -rh ${params.ldpop_rho_range} --approx > lookupTable.txt
    """
}


process LDHAT_CONVERT{
//     publishDir "Output", mode: "copy", saveAs: {filename -> "${fn_modifier}_${filename}"}

    input:
        tuple val(mutation_rate),
            val(sample_size),
            path("LDhat_reformated.fa"),
            path("lookupTable.txt")

    output:
        tuple val(mutation_rate),
            val(sample_size),
            path("lookupTable.txt"),
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


process SWITCH_TO_GENE_CONVERSION_MODE{
    publishDir "Output", mode: "copy", saveAs: {filename -> "${fn_modifier}_${filename}"}

    input:
        tuple val(mutation_rate),
            val(sample_size),
            path("lookupTable.txt"),
            path("locs.txt"),
            path("sites.txt")

    output:
        tuple val(mutation_rate),
            val(sample_size),
            path("lookupTable.txt"),
            path("sites.txt"),
            path("locs_C.txt")

    script:
        """
        get_ldhat_to_use_gene_conv.py locs.txt
        """

}

process LDHAT_PAIRWISE{
    publishDir "Output", mode: "copy", saveAs: {filename -> "${fn_modifier}_${filename}"}

    input:
        tuple val(mutation_rate),
            val(sample_size),
            path("lookupTable.txt"),
            path("sites.txt"),
            path("locs_C.txt")

    output:
        tuple val(mutation_rate),
            val(sample_size),
            path("pairwise_freqs.txt"),
            path("pairwise_outfile.txt"),
            path("pairwise_stdOut.txt")

    script:
        // uses pexpect to handle unavoidable prompts
        """
        run_pairwise_with_pexpect.py ${params.recom_tract_len} sites.txt locs_C.txt lookupTable.txt > pairwise_stdOut.txt
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
    params.recom_tract_len = 500
    params.ldpop_rho_range = "101,100"

    params.input_fasta = 'none'

    // Input verification
    if (params.input_fasta == 'none') {
        println "No input .fa specified. Use --input_fasta [.fa]"
        exit 1
    }

    // Channels
    input_fasta_channel = Channel.fromPath( params.input_fasta )

    // For each process there is a output of tuple with the params that change + necessary files/values  to move forward until they are no longer need

    LDHAT_REFORMAT_FASTA(input_fasta_channel, params.mutation_rate)

    LOOKUP_TABLE_LDPOP(LDHAT_REFORMAT_FASTA.out)

    LDHAT_CONVERT(LOOKUP_TABLE_LDPOP.out)

    SWITCH_TO_GENE_CONVERSION_MODE(LDHAT_CONVERT.out)

    LDHAT_PAIRWISE(SWITCH_TO_GENE_CONVERSION_MODE.out)

    // PAIRWISE_PROCESS_OUTPUT(LDHAT_PAIRWISE.out)

    // collectedFile = PAIRWISE_PROCESS_OUTPUT.out.processed_results_csv.collectFile(name:"collected_results.csv",storeDir:"Output/Results", keepHeader:true)

    // PLOT_RESULTS(collectedFile)
}