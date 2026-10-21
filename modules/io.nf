/*
 * Input stages shared by both arms. Copied verbatim from main.nf (process names and scripts
 * unchanged, so the work cache is shared and a faithful run's downloads are reused).
 */
nextflow.enable.dsl = 2

process ENA_SAMPLESHEET {
    // Ask the European Nucleotide Archive which FASTQ files belong to the project.
    publishDir "${params.outdir}/00_samplesheet", mode: 'copy'
    conda 'conda-forge::python=3.11'
    input:  val accession
    output: path 'samplesheet.csv'
    script: "ena_samplesheet.py ${accession} -o samplesheet.csv --include '${params.samples}'"
    stub:   "printf 'sample,run,seq_run,fastq_1,fastq_2,md5_1,md5_2,read_count\\nS1,R1,run_1,u1,u2,m1,m2,10\\nS2,R2,run_2,u1,u2,m1,m2,10\\n' > samplesheet.csv"
}

process DOWNLOAD {
    // Fetch one paired-end run and verify the checksums published by the archive.
    tag "$sample"
    maxForks 4
    errorStrategy 'retry'; maxRetries 3
    input:  tuple val(sample), val(url1), val(url2), val(md5_1), val(md5_2)
    output: tuple val(sample), path("${sample}_1.fastq.gz"), path("${sample}_2.fastq.gz")
    script:
    """
    curl -sSL --retry 5 -o ${sample}_1.fastq.gz ${url1}
    curl -sSL --retry 5 -o ${sample}_2.fastq.gz ${url2}
    echo "${md5_1}  ${sample}_1.fastq.gz" | md5sum -c -
    echo "${md5_2}  ${sample}_2.fastq.gz" | md5sum -c -
    """
    stub: "touch ${sample}_1.fastq.gz ${sample}_2.fastq.gz"
}
