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
            path("theta_est.csv")

    script:
    """
    fasta_variant_sites.py LDhat_reformated.fa
    snps=\$(wc -l variants_in_fasta.csv | awk '{split(\$0,a," "); print a[1]}')
    genome_len=\$(head -1 locs_C.txt | awk '{split(\$0,a," "); print a[2]}' | awk '{split(\$0,a,"."); print a[1]}')
    original_watterson.py \$genome_len \$snps ${sample_size} > theta_est.csv
    """
}



workflow {
    // Note: Channels can be called unlimited number of times in DSL2
    // A process component can be invoked only once in the same workflow context

    params.prefix_filename = 'none'
    params.input_fasta = 'none'


    // Input verification
    if (params.input_fasta == 'none') {
        println "No input .fa specified. Use --input_fasta [.fa]"
        exit 1
    }

    // Channels
    input_fasta_channel = Channel.fromPath( params.input_fasta )

    // For each process there is a output of tuple with the params that change + necessary files/values  to move forward until they are no longer need
    PREFIX_FILENAME(input_fasta_channel, params.prefix_filename)

    LDHAT_REFORMAT_FASTA(PREFIX_FILENAME.out)

    LDHAT_CONVERT(LDHAT_REFORMAT_FASTA.out)

    SWITCH_TO_GENE_CONVERSION_MODE(LDHAT_CONVERT.out)

    WATTERSON_ESTIMATE(SWITCH_TO_GENE_CONVERSION_MODE.out)
}