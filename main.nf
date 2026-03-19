nextflow.enable.dsl=2

params.reads = "$projectDir/data/raw/*_{1,2}.fastq.gz"
params.outdir = "$projectDir/results"

process FASTQC {
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "*_fastqc.{zip,html}"

    script:
    """
    fastqc -q $reads
    """
}

workflow {
    read_pairs_ch = channel.fromFilePairs(params.reads, checkIfExists: true)
    FASTQC(read_pairs_ch)
}