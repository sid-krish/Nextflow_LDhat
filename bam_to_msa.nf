#! /usr/bin/env nextflow
nextflow.enable.dsl = 2


process PREPEND_FILENAME {

    maxForks 1

    // echo true

    input:
        tuple path(bam),
            path(fasta)
        
        val(prepend_fn)

    output:
        tuple stdout,
            path(bam),
            path(fasta)

    script:
    """
    prepend_filename_ldhat.py ${bam} ${prepend_fn} 
    """
}


process FREEBAYES {
    publishDir "Bam_to_msa_Output", mode: "copy", saveAs: {filename -> "${prepend_filename}${filename}"}

    input:
        tuple val(prepend_filename),
            path(bam),
            path(fasta)

    output:
        tuple val(prepend_filename),
            path("freeBayesOut.vcf"),
            path(fasta)

    script:
    """
    samtools faidx ${fasta}
    samtools sort --threads $task.cpus ${bam} -o Aligned.csorted.bam
    samtools index -@ $task.cpus Aligned.csorted.bam

    freebayes -f ${fasta} -p 1 Aligned.csorted.bam > freeBayesOut.vcf
    """
}


process VCF_TO_FASTA {
    input:

    output:

    script:
    """

    """
}


process UPDATE_HEADER {
    input:

    output:

    script:
    """

    """
}


process MERGE_FASTAS {
    input:

    output:

    script:
    """

    """
}


workflow {
    // Params
    params.prepend_filename = "none"
    params.bam_file = 'none'
    params.reference_genome = 'none'

    // Channels
    bam_file_channel = Channel.fromPath( params.bam_file )
    reference_genome_channel = Channel.fromPath( params.reference_genome )

    bam_and_fa = bam_file_channel.combine(reference_genome_channel)
    

    // Process execution
    PREPEND_FILENAME(bam_and_fa, params.prepend_filename)

    FREEBAYES(PREPEND_FILENAME.out)

}