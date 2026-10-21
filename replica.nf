#!/usr/bin/env nextflow
/*
 * replica: a portable re-implementation of the NBGW plant metabarcoding pipeline
 * (github.com/colford/nbgw-plant-illumina-pipeline), same tools, parameters and order.
 * One process per stage, in the order the reads flow. See docs/replica.md.
 */
nextflow.enable.dsl = 2

include { ENA_SAMPLESHEET; DOWNLOAD } from './modules/io.nf'
include { REFERENCE } from './modules/reference.nf'

// 0. the reads come from modules/io.nf
// 1. trim the reads
process TRIM {
    // Trimmomatic, their exact recipe: adapters, HEADCROP:3, 4-base sliding window, MINLEN 200.
    // Only intact pairs go on.
    tag "$sample"
    conda params.conda_trimmomatic
    input:  tuple val(sample), path(r1), path(r2)
    output: tuple val(sample), path("${sample}.P1.fastq.gz"), path("${sample}.P2.fastq.gz"), emit: reads
            path "${sample}.trim_stats.tsv", emit: stats
    script:
    """
    # a path of your own, else the stock list of that name inside the installed Trimmomatic
    ad='${params.adapters}'
    [ -f "\$ad" ] || ad=\$(ls \$(dirname \$(command -v trimmomatic))/../share/trimmomatic*/adapters/"\$ad" | head -1)
    [ -f "\$ad" ] || { echo "no adapter list '${params.adapters}'" >&2; exit 1; }

    trimmomatic PE -threads ${task.cpus} -phred33 ${r1} ${r2} \\
        ${sample}.P1.fastq.gz ${sample}.U1.fastq.gz ${sample}.P2.fastq.gz ${sample}.U2.fastq.gz \\
        ILLUMINACLIP:\$ad:2:30:10 HEADCROP:3 LEADING:20 \\
        SLIDINGWINDOW:4:${params.trim_quality} MINLEN:200
    raw=\$(( \$(zcat ${r1} | wc -l) / 4 )); kept=\$(( \$(zcat ${sample}.P1.fastq.gz | wc -l) / 4 ))
    printf 'sample\\traw_pairs\\ttrimmed_pairs\\n${sample}\\t%s\\t%s\\n' \$raw \$kept > ${sample}.trim_stats.tsv
    """
    stub: "touch ${sample}.P1.fastq.gz ${sample}.P2.fastq.gz; printf 'sample\\traw_pairs\\ttrimmed_pairs\\n${sample}\\t10\\t9\\n' > ${sample}.trim_stats.tsv"
}

// 2. merge the pairs
process MERGE {
    // FLASH joins each forward/reverse pair into one sequence spanning the whole amplicon.
    tag "$sample"
    conda params.conda_flash
    input:  tuple val(sample), path(p1), path(p2)
    output: tuple val(sample), path("${sample}.extendedFrags.fastq.gz"), emit: merged
    script:
    """
    flash --max-overlap 450 --min-overlap 10 --max-mismatch-density 0.25 \\
          -t ${task.cpus} -z -o ${sample} ${p1} ${p2}
    """
    stub: "touch ${sample}.extendedFrags.fastq.gz"
}

// 3. split by marker and length
process SPLIT_BY_PRIMER {
    // The two markers were pooled before sequencing, so a merged read is either rbcL or ITS2:
    // told apart by a primer fragment at each end (primers stay on), then the minimum merged
    // length of their script (450 rbcL, 300 ITS2).
    tag "$sample"
    conda 'conda-forge::python=3.11'
    input:  tuple val(sample), path(merged)
    output: tuple val(sample), path("${sample}.*.fasta"), emit: fasta
            path "${sample}.split_stats.tsv", emit: stats
    script: "split_by_primer.py --sample ${sample} --merged ${merged} --mismatches ${params.primer_mismatches} --markers '${groovy.json.JsonOutput.toJson(params.markers)}'"
    stub:   "for m in ${params.markers.keySet().join(' ')}; do touch ${sample}.\$m.fasta; done; printf 'sample\\tmerged\\n${sample}\\t8\\n' > ${sample}.split_stats.tsv"
}

// 4. dereplicate within each sample
process DEREP {
    // Identical reads collapse to one with ;size=N. Labels become SAMPLE|n so the sample
    // survives the pooling in the next step.
    tag "$sample:$marker"
    conda params.conda_vsearch
    input:  tuple val(sample), val(marker), path(fasta)
    output: tuple val(marker), path("${sample}.${marker}.derep.fasta")
    script:
    """
    if [ -s ${fasta} ]; then
      vsearch --derep_fulllength ${fasta} --sizeout --relabel '${sample}|' --fasta_width 0 --output ${sample}.${marker}.derep.fasta
    else touch ${sample}.${marker}.derep.fasta; fi
    """
    stub: "touch ${sample}.${marker}.derep.fasta"
}

// 5. cluster across samples, drop singletons
process CLUSTER {
    // All samples pooled and clustered at 100% identity. A sequence seen once in the whole
    // study is treated as an error and dropped.
    tag "$marker"
    publishDir "${params.outdir}/05_clusters", mode: 'copy'
    conda params.conda_vsearch
    input:  tuple val(marker), path(fastas)
    output: tuple val(marker), path("${marker}.centroids.fasta"), path("${marker}.uc")
    script:
    """
    cat ${fastas} > all.fasta
    if [ "${params.cluster_mode}" = "align" ]; then
      # their exact command: 100% identity by alignment (terminal gaps ignored). Slow: O(n^2).
      vsearch --cluster_size all.fasta --id 1.0 --sizein --sizeout --fasta_width 0 \\
              --uc ${marker}.uc --centroids all.centroids.fasta --threads ${task.cpus}
    else
      # default: identical sequences pooled across samples. Same clusters except that a sequence
      # contained inside a longer one is kept apart instead of absorbed. Seconds instead of hours.
      vsearch --derep_fulllength all.fasta --sizein --sizeout --fasta_width 0 \\
              --uc ${marker}.uc --output all.centroids.fasta
    fi
    vsearch --sortbysize all.centroids.fasta --minsize 2 --fasta_width 0 --output ${marker}.centroids.fasta
    """
    stub: "printf '>c1\\nACGT\\n>c2\\nACGT\\n' > ${marker}.centroids.fasta; touch ${marker}.uc"
}

// 6. the reference library is in modules/reference.nf
// 7. blast the clusters and assign a taxon
process BLAST {
    // Top 20 hits per cluster. Chunked so an interrupted run loses at most one chunk; primers
    // are trimmed off the query first (trim_primers.py).
    tag "$marker:${chunk.baseName}"
    cache 'deep'
    conda "conda-forge::python=3.11 ${params.conda_blast}"
    input:  tuple val(marker), path(chunk)
            path refdb
    output: tuple val(marker), path("${marker}.${chunk.baseName}.blast.tsv")
    script:
    """
    trim_primers.py --fasta ${chunk} --marker ${marker} --mismatches ${params.primer_mismatches} \\
        --markers '${groovy.json.JsonOutput.toJson(params.markers)}' > query.fasta 2> trim.log
    blastn -query query.fasta -db ${refdb}/${marker} -max_target_seqs 20 -outfmt '6 std' \\
           -num_threads ${task.cpus} -out ${marker}.${chunk.baseName}.blast.tsv
    """
    stub: "touch ${marker}.${chunk.baseName}.blast.tsv"
}

process ASSIGN {
    // The paper's rules: agreement among the best-scoring hits decides species, genus or family.
    // The 1st-percentile bit-score cut-off is computed over all clusters of the marker at once.
    tag "$marker"
    publishDir "${params.outdir}/07_assignment", mode: 'copy', pattern: '*.assignments.tsv'
    conda 'conda-forge::python=3.11'
    input:  tuple val(marker), path(blasts), path(uc)
            path refdb
    output: tuple val(marker), path(uc), path("${marker}.assignments.tsv")
    script:
    """
    cat ${blasts} > ${marker}.blast.tsv
    assign_taxonomy.py --blast ${marker}.blast.tsv --taxonomy ${refdb}/${marker}.taxonomy.tsv --out ${marker}.assignments.tsv
    """
    stub: "printf 'cluster\\ttaxon\\trank\\ttop_bitscore\\ttop_identity\\tn_top_hits\\tnote\\n' > ${marker}.assignments.tsv"
}

// 8. the tables
process TABLES {
    // Sum reads per cluster and per taxon for every sample: the sample-by-taxon matrix.
    tag "$marker"
    publishDir "${params.outdir}/08_tables", mode: 'copy'
    conda 'conda-forge::python=3.11'
    input:  tuple val(marker), path(uc), path(assignments)
    output: tuple val(marker), path("${marker}.taxon_table.tsv"), emit: taxa
            path "${marker}.cluster_table.tsv"
            path "${marker}.table_stats.tsv", emit: stats
    script: "build_tables.py --uc ${uc} --assignments ${assignments} --marker ${marker}"
    stub:   "touch ${marker}.taxon_table.tsv ${marker}.cluster_table.tsv; printf 'sample\\t${marker}_assigned\\nS1\\t5\\n' > ${marker}.table_stats.tsv"
}

// 9. compare with the publication
process COMPARE {
    // Both ranks, because a label-rank disagreement and a family-rank one mean different things:
    // the first can be two names for one plant, the second cannot.
    tag "$marker:$rank"
    publishDir "${params.outdir}/09_comparison", mode: 'copy'
    conda 'conda-forge::python=3.11 conda-forge::openpyxl'
    input:  tuple val(marker), path(table)
            path published
            path refdb
            each rank
    output: path "${marker}.vs_published.${rank}*.tsv"
    script: """
    compare_published.py --replica ${table} --published ${published} --sheet ${params.published_sheets[marker]} \\
        --marker ${marker} --taxonomy ${refdb}/${marker}.taxonomy.tsv --rank ${rank} \\
        --out ${marker}.vs_published.${rank}.tsv
    """
    stub:   "touch ${marker}.vs_published.${rank}.tsv ${marker}.vs_published.${rank}.name_translation.tsv"
}

process READ_TRACKING {
    // One row per sample, one column per stage: where reads are lost, and whether totals match the paper.
    publishDir "${params.outdir}", mode: 'copy'
    conda 'conda-forge::python=3.11'
    input:  path stats
    output: path 'read_tracking.tsv'
    script: "merge_stats.py ${stats} > read_tracking.tsv"
    stub:   "touch read_tracking.tsv"
}

workflow {
    sheet = params.input ? Channel.fromPath(params.input) : ENA_SAMPLESHEET(params.accession)
    def limit = (params.max_samples as Integer)
    runs  = sheet.splitCsv(header: true)
                 .take(limit)
                 .map { r -> tuple(r.sample, r.fastq_1, r.fastq_2, r.md5_1, r.md5_2) }

    DOWNLOAD(runs)
    TRIM(DOWNLOAD.out)
    MERGE(TRIM.out.reads)
    SPLIT_BY_PRIMER(MERGE.out.merged)

    per_marker = SPLIT_BY_PRIMER.out.fasta
        .flatMap { sample, files -> [files].flatten().collect { f -> tuple(sample, f.name.tokenize('.')[-2], f) } }
    DEREP(per_marker)
    CLUSTER(DEREP.out.groupTuple(sort: { it.name }))

    REFERENCE(file(params.ref_manifest, checkIfExists: true),
              file(params.ref_sequences, checkIfExists: true))
    chunks = CLUSTER.out.map { m, cent, uc -> tuple(m, cent) }
                        .splitFasta(by: (params.blast_chunk as Integer).intValue(), file: true, elem: 1)
    BLAST(chunks, REFERENCE.out.db)
    hits = BLAST.out.groupTuple()
                    .join(CLUSTER.out.map { m, cent, uc -> tuple(m, uc) })
    ASSIGN(hits, REFERENCE.out.db)
    TABLES(ASSIGN.out)

    if (params.published)
        COMPARE(TABLES.out.taxa, file(params.published, checkIfExists: true), REFERENCE.out.db,
                ['label', 'family'])

    READ_TRACKING(TRIM.out.stats.mix(SPLIT_BY_PRIMER.out.stats, TABLES.out.stats).collect())
}
