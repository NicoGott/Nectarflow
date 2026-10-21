#!/usr/bin/env nextflow
/*
 * replica, second arm: the same study analysed with today's tools. replica.nf stays the faithful
 * re-implementation of the authors' pipeline; this arm is where the method differs. Six phases:
 * cutadapt, DADA2, classify, controls, compare (with the paper and with arm 1), ecology.
 * See docs/reanalysis.md.
 */
nextflow.enable.dsl = 2

include { ENA_SAMPLESHEET; DOWNLOAD } from './modules/io.nf'
include { REFERENCE } from './modules/reference.nf'

// 1. cut the spacer and the primers, split the markers
process CUTADAPT {
    // One linked adapter per read: 5' anchored and required (6N spacer + primer), 3' optional
    // (the other primer, reached only when the amplicon is shorter than the read, as for ITS2).
    // --discard-untrimmed with --pair-filter=any makes trimming and marker assignment one
    // decision. Arm 1 leaves the spacer and the primers on; here the denoiser sees neither.
    tag "$sample:$marker"
    conda params.conda_cutadapt
    input:  tuple val(sample), path(r1), path(r2)
            each marker
    output: tuple val(sample), val(marker), path("${sample}.${marker}.R1.fastq.gz"), path("${sample}.${marker}.R2.fastq.gz"), emit: reads
            path "${sample}.${marker}.cutadapt.json", emit: stats
    script:
    def m = params.primers[marker]
    """
    cutadapt -j ${task.cpus} -e ${params.primer_mismatches} --no-indels \\
        -a '^${m.fwd}...${m.rc_rev}' \\
        -A '^${m.rev}...${m.rc_fwd}' \\
        --discard-untrimmed --pair-filter=any --minimum-length ${params.min_len_after_trim} \\
        --json ${sample}.${marker}.cutadapt.json \\
        -o ${sample}.${marker}.R1.fastq.gz -p ${sample}.${marker}.R2.fastq.gz \\
        ${r1} ${r2}
    """
    stub: "touch ${sample}.${marker}.R1.fastq.gz ${sample}.${marker}.R2.fastq.gz; echo '{}' > ${sample}.${marker}.cutadapt.json"
}

process CUTADAPT_STATS {
    // One table: pairs in, pairs kept per marker, and what neither primer pair claimed.
    publishDir "${params.outdir}/01_cutadapt", mode: 'copy'
    conda 'conda-forge::python=3.11'
    input:  path jsons
    output: path 'cutadapt_stats.tsv'
    script: "cutadapt_stats.py --out cutadapt_stats.tsv ${jsons}"
    stub:   "touch cutadapt_stats.tsv"
}

// 2. denoise into ASVs
process DADA2_DENOISE {
    // One task per marker per sequencing run: the error model is a property of the run (seq_run
    // comes from the archive's library_name, run_1 = 2018, run_2 and run_3 = 2019). Filtering
    // happens here so DADA2 sees the original quality profiles. rbcL is truncated (fixed 458 bp
    // insert, R2 quality falls after ~210), ITS2 is not, its length varies between plants.
    // pool='pseudo' recovers plants present in a single insect.
    tag "$marker:$run"
    conda params.conda_dada2
    cpus 4
    input:  tuple val(marker), val(run), path(reads)
    output: tuple val(marker), path("${marker}.${run}.seqtab.rds"),      emit: seqtab
            tuple val(marker), path("${marker}.${run}.dada2_track.tsv"), emit: track
    script:
    def m = params.dada2[marker]
    """
    dada2_denoise.R --marker=${marker} --run=${run} --trunc_len=${m.trunc_len} --max_ee=${m.max_ee} \\
        --min_overlap=${m.min_overlap} --threads=${task.cpus} ${reads}
    """
    stub: "touch ${marker}.${run}.seqtab.rds ${marker}.${run}.dada2_track.tsv"
}

process DADA2_CHIMERA {
    // Chimeras are made in the PCR, so they are looked for across the whole study, not within a run.
    tag "$marker"
    publishDir "${params.outdir}/02_asv", mode: 'copy'
    conda params.conda_dada2
    cpus 4
    input:  tuple val(marker), path(tabs), path(tracks)
    output: tuple val(marker), path("${marker}.asv_table.tsv"), path("${marker}.asvs.fasta"), emit: asv
            path "${marker}.dada2_track.tsv", emit: track
    script: "dada2_chimera.R --marker=${marker} --threads=${task.cpus} ${tabs} ${tracks}"
    stub:   "touch ${marker}.asv_table.tsv ${marker}.asvs.fasta ${marker}.dada2_track.tsv"
}

process CURATE {
    // Post-clustering curation, the LULU method of Froslev et al. 2017 (Nat Commun 8:1188), which
    // they validated on plant ITS2 exactly like this. The rDNA is multi-copy and its copies inside
    // one genome differ, so one plant yields several ASVs a base or two apart; DADA2 keeps them,
    // rightly, because they are real sequences, but they split one plant across several rows.
    // An ASV found almost only where a more widespread, similar ASV is found, and always rarer
    // there, is merged into it. No reference library is involved, so this cannot invent agreement.
    // Run with mumu, Mahe's maintained C++ implementation of the method, at its defaults, which are
    // the published ones. mumu documents three deliberate differences from the 2017 R package, the
    // largest being that lulu's handling of a partial overlap cancelled its own co-occurrence
    // option (lulu issue 8); --legacy reproduces the R behaviour if that is ever wanted.
    tag "$marker"
    publishDir "${params.outdir}/02_asv", mode: 'copy'
    conda "conda-forge::python=3.11 ${params.conda_vsearch_current} ${params.conda_mumu}"
    cpus 4
    input:  tuple val(marker), path(table), path(fasta)
    output: tuple val(marker), path("${marker}.lulu.asv_table.tsv"), path("${marker}.lulu.asvs.fasta"), emit: asv
            path "${marker}.lulu_merges.tsv"
            path "${marker}.mumu_log.tsv"
    script:
    def l = params.lulu
    """
    vsearch --usearch_global ${fasta} --db ${fasta} --self --id ${l.min_match / 100} --iddef 1 \\
        --userout match_list.txt --userfields query+target+id --maxaccepts 0 --query_cov .9 \\
        --maxhits 10 --threads ${task.cpus}

    # mumu reads OTUs in rows; this pipeline carries samples in rows, so turn it over and back
    transpose.py < ${table} > otus_in_rows.tsv
    mumu --otu_table otus_in_rows.tsv --match_list match_list.txt \\
        --minimum_match ${l.min_match} --minimum_relative_cooccurrence ${l.min_cooccurrence} \\
        --minimum_ratio ${l.min_ratio} --minimum_ratio_type min \\
        --new_otu_table curated_otus_in_rows.tsv --log ${marker}.mumu_log.tsv
    transpose.py < curated_otus_in_rows.tsv > ${marker}.lulu.asv_table.tsv

    # the ASVs that survived, and their sequences, for the BLAST that follows
    tail -n +2 curated_otus_in_rows.tsv | cut -f1 > kept.txt
    vsearch --fastx_getseqs ${fasta} --labels kept.txt --fastaout ${marker}.lulu.asvs.fasta

    # every merge mumu made, so each one can be checked: columns are query, parent, the query's
    # total reads and the number of samples it was in (log columns 1, 2, 4 and 8)
    awk -F'\\t' 'BEGIN{OFS="\\t"; print "asv","merged_into","reads","samples"} \\
        \$18=="accepted"{print \$1,\$2,\$4,\$8}' ${marker}.mumu_log.tsv > ${marker}.lulu_merges.tsv
    """
    stub: "touch ${marker}.lulu.asv_table.tsv ${marker}.lulu.asvs.fasta ${marker}.lulu_merges.tsv ${marker}.mumu_log.tsv"
}

// 3. classify the ASVs
process BLAST_ASV {
    // Same reference library as arm 1, so a difference between the arms is one of method. No
    // primer trimming here: phase 1 removed the primers, so a hit cannot be scored on them, as it
    // could in arm 1 where GenBank references carrying the primer outscored Barcode UK ones.
    tag "$marker"
    conda params.conda_blast_current
    cpus 4
    input:  tuple val(marker), path(table), path(fasta)
            path refdb
    output: tuple val(marker), path(table), path(fasta), path("${marker}.blast.tsv"), path(refdb)
    script:
    """
    blastn -query ${fasta} -db ${refdb}/${marker} -max_target_seqs ${params.max_target_seqs} \\
           -outfmt '6 qseqid sseqid pident length qlen qcovhsp mismatch gapopen evalue bitscore score' \\
           -num_threads ${task.cpus} -out ${marker}.blast.tsv
    """
    stub: "touch ${marker}.blast.tsv"
}

process ASSIGN_ASV {
    // Rank floors calibrated on this library, tie band in substitutions, lowest common ancestor
    // over the candidates; a rule may only make a call coarser. See bin/assign_taxonomy_asv.py.
    tag "$marker"
    publishDir "${params.outdir}/03_assignment", mode: 'copy'
    conda 'conda-forge::python=3.11'
    input:  tuple val(marker), path(table), path(fasta), path(blast), path(refdb)
    output: tuple val(marker), path("${marker}.taxon_table.tsv"), emit: taxa
            path "${marker}.asv_taxonomy.tsv", emit: assignments
            path "${marker}.taxon_table.full.tsv"
            path "${marker}.rank_summary.tsv"
            path "${marker}.classify_params.tsv"
    script:
    """
    assign_taxonomy_asv.py --blast ${blast} --taxonomy ${refdb}/${marker}.taxonomy.tsv \\
        --asv-table ${table} --asv-fasta ${fasta} --marker ${marker} \\
        --max-target-seqs ${params.max_target_seqs} --out .
    """
    stub: "touch ${marker}.taxon_table.tsv ${marker}.asv_taxonomy.tsv ${marker}.taxon_table.full.tsv ${marker}.rank_summary.tsv ${marker}.classify_params.tsv"
}

process CONSENSUS {
    // The two markers combined into one composition per insect: the operation the authors did by
    // hand, but from our own reference lineages and never from their mapping sheet. Arm 1 never
    // combined its markers either, and phase 6 needs one table per arm, so the same process runs
    // on both arms' marker tables and publishes under the arm's name.
    tag "$arm"
    publishDir { "${params.outdir}/04_consensus/${arm}" }, mode: 'copy'
    conda 'conda-forge::python=3.11'
    input:  tuple val(arm), path(rbcl), path(its2)
            path refdb
    output: tuple val(arm), path("consensus_taxon_table.tsv"), emit: table
            path "consensus_*.tsv"
    script:
    """
    consensus_markers.py --rbcl ${rbcl} --its2 ${its2} \\
        --taxonomy ${refdb}/rbcL.taxonomy.tsv ${refdb}/ITS2.taxonomy.tsv --out .
    """
    stub: "touch consensus_taxon_table.tsv consensus_map.tsv"
}

// 4. what is in the blanks
process CONTROLS {
    // Reports, never subtracts: a blank holds contamination and index bleed, and only a reader
    // can tell which one a given ASV is.
    publishDir "${params.outdir}/05_controls", mode: 'copy'
    conda 'conda-forge::python=3.11'
    input:  path tables
            path assignments
    output: path "controls_*"
    script:
    """
    controls_report.py --asv-table ${tables} --assignments ${assignments} \\
        --controls '${params.controls}' --out .
    """
    stub: "touch controls_summary.txt"
}

// 5. compare, with the publication and with arm 1
process COMPARE_PUBLISHED_ARM2 {
    // The same reader as arm 1 (bin/compare_published.py), so both arms' agreement with the
    // Dryad sheets comes from one code path.
    tag "$marker:$rank"
    publishDir "${params.outdir}/06_comparison", mode: 'copy'
    conda 'conda-forge::python=3.11 conda-forge::openpyxl'
    input:  tuple val(marker), path(table)
            path published
            path refdb
            each rank
    output: tuple val(marker), val(rank), path("${marker}.arm2_vs_published.${rank}.tsv"),
                  path("${marker}.arm2_vs_published.${rank}.name_translation.tsv"), emit: vs_pub
            path "${marker}.arm2_vs_published.${rank}*.tsv"
    script:
    """
    compare_published.py --replica ${table} --published ${published} --sheet ${params.published_sheets[marker]} \\
        --marker ${marker} --taxonomy ${refdb}/${marker}.taxonomy.tsv --rank ${rank} \\
        --out ${marker}.arm2_vs_published.${rank}.tsv
    """
    stub: "touch ${marker}.arm2_vs_published.${rank}.tsv ${marker}.arm2_vs_published.${rank}.name_translation.tsv"
}

process COMPARE_ARMS {
    // Arm 1 against arm 2 per insect, at label level and rolled up to family.
    tag "$marker:$rank"
    publishDir "${params.outdir}/06_comparison", mode: 'copy'
    conda 'conda-forge::python=3.11'
    // At the label rank --translate1/2 put both arms in the authors' vocabulary, so this row reads
    // next to the two published rows; without them a coarser arm scores low for using its own names.
    // At the family rank compare_arms.py rolls both arms up instead and ignores them.
    // Arm 1's three files are declared inputs, not interpolated paths, so they are staged and hashed.
    input:  tuple val(marker), path(arm2_table), val(rank), path(vs_pub), path(vs_pub_translation),
                  path(arm1, stageAs: 'arm1/*'), path(arm1_pub, stageAs: 'arm1/*'),
                  path(replica_translation, stageAs: 'arm1/*')
            path refdb
    output: path "${marker}.arms.${rank}*.tsv"
    script:
    """
    compare_arms.py --arm1 ${arm1} --arm2 ${arm2_table} --marker ${marker} --rank ${rank} \\
        --translate1 ${replica_translation} --translate2 ${vs_pub_translation} \\
        --lineage ${refdb}/rbcL.taxonomy.tsv ${refdb}/ITS2.taxonomy.tsv \\
        --vs-published-arm1 ${arm1_pub} --vs-published-arm2 ${vs_pub} \\
        --out ${marker}.arms.${rank}.tsv --summary ${marker}.arms.${rank}.summary.tsv
    """
    stub: "touch ${marker}.arms.${rank}.tsv ${marker}.arms.${rank}.summary.tsv"
}

// 6. the ecology
process ECOLOGY_METADATA {
    // The one fact a sample code cannot carry: which insect the pollen came off. From the
    // metadata sheet only, never the published read counts that phase 5 compares against.
    publishDir "${params.outdir}/07_ecology", mode: 'copy'
    conda 'conda-forge::python=3.11 conda-forge::openpyxl'
    input:  path published
            path samplesheet        // carries seq_run, so a study variable that is a function of
                                    // the sequencing batch can be reported as one
    output: path 'insect_metadata.tsv'
    script:
    def sheet = samplesheet ? "--samplesheet ${samplesheet} " : ""
    "ecology_metadata.py --published ${published} ${sheet}--out insect_metadata.tsv"
    stub:   "printf 'sample\\tinsect_group\\n' > insect_metadata.tsv"
}

process PUBLISHED_TABLE {
    // The authors' combined table in our shape, read counts copied and nothing else.
    conda 'conda-forge::python=3.11 conda-forge::openpyxl'
    input:  path published
    output: path 'published_taxon_table.tsv'
    script: "published_table.py --published ${published} --out published_taxon_table.tsv"
    stub:   "printf 'taxon\n' > published_taxon_table.tsv"
}

process ECOLOGY {
    // One script, three tables (arm 1, arm 2, published), so the arms are compared on their
    // conclusions and not only on their read counts.
    tag "$label"
    publishDir { "${params.outdir}/07_ecology/${label}" }, mode: 'copy'
    conda params.conda_ecology
    // ecology.R is single-threaded: mvabund's anova has no parallel argument and vegan is only
    // asked for one adonis2. Declaring more would hold cores it cannot use, which on a 4-core box
    // froze the whole denoising queue behind one arm for hours.
    cpus 1
    input:  tuple val(label), path(table)
            path metadata
    output: path '*.tsv'
            path 'summary.txt'
    script:
    """
    ecology.R --table=${table} --metadata=${metadata} --label=${label} --out=.
    """
    stub: "touch summary.txt run_info.tsv"
}

workflow {
    // .first() makes this a value channel: the samplesheet is read three times below, and a queue
    // channel would be emptied by the first reader, which left ECOLOGY_METADATA without a seq_run
    // column and so with no way to notice that year and sequencing run are the same thing here.
    sheet = (params.input ? Channel.fromPath(params.input) : ENA_SAMPLESHEET(params.accession)).first()
    def limit = (params.max_samples as Integer)
    runs  = sheet.splitCsv(header: true)
                 .take(limit)
                 .map { r -> tuple(r.sample, r.fastq_1, r.fastq_2, r.md5_1, r.md5_2) }

    seq_run = sheet.splitCsv(header: true)
                 .map { r -> tuple(r.sample, r.seq_run ?: 'run_1') }

    DOWNLOAD(runs)
    CUTADAPT(DOWNLOAD.out, params.primers.keySet() as List)
    CUTADAPT_STATS(CUTADAPT.out.stats.collect())

    // group the trimmed pairs by (marker, sequencing run)
    by_run = CUTADAPT.out.reads
        .combine(seq_run, by: 0)
        .map { sample, marker, r1, r2, run -> tuple(tuple(marker, run), [r1, r2]) }
        .groupTuple()
        .map { key, files -> tuple(key[0], key[1], files.flatten()) }

    DADA2_DENOISE(by_run)
    DADA2_CHIMERA(DADA2_DENOISE.out.seqtab.groupTuple().join(DADA2_DENOISE.out.track.groupTuple()))

    REFERENCE(file(params.ref_manifest, checkIfExists: true),
              file(params.ref_sequences, checkIfExists: true))
    CURATE(DADA2_CHIMERA.out.asv)
    BLAST_ASV(CURATE.out.asv, REFERENCE.out.db)
    ASSIGN_ASV(BLAST_ASV.out)

    taxa = ASSIGN_ASV.out.taxa.branch { marker, t -> rbcL: marker == 'rbcL'; ITS2: marker == 'ITS2' }
    // one consensus process for both arms: arm 2's own marker tables, and arm 1's if a faithful
    // run has been pointed at with --replica_results
    markers = taxa.rbcL.map { m, t -> t }.combine(taxa.ITS2.map { m, t -> t })
                  .map { r, i -> tuple('arm2', r, i) }
    if (params.replica_results) {
        markers = markers.mix(Channel.of(tuple('arm1',
            file("${params.replica_results}/08_tables/rbcL.taxon_table.tsv", checkIfExists: true),
            file("${params.replica_results}/08_tables/ITS2.taxon_table.tsv", checkIfExists: true))))
    }
    CONSENSUS(markers, REFERENCE.out.db)

    CONTROLS(CURATE.out.asv.map { m, t, f -> t }.collect(),
             ASSIGN_ASV.out.assignments.collect())

    if (params.published) {
        COMPARE_PUBLISHED_ARM2(ASSIGN_ASV.out.taxa, file(params.published, checkIfExists: true),
                               REFERENCE.out.db, ['label', 'family'])
        if (params.replica_results) {
            // each rank of the arms comparison gets the published comparison made at the same rank,
            // with the name translation that run wrote, and arm 1's three files alongside them
            pub = COMPARE_PUBLISHED_ARM2.out.vs_pub
                .map { marker, rank, f, translation -> tuple(marker, rank, f, translation,
                       file("${params.replica_results}/08_tables/${marker}.taxon_table.tsv", checkIfExists: true),
                       file("${params.replica_results}/09_comparison/${marker}.vs_published.${rank}.tsv", checkIfExists: true),
                       file("${params.replica_results}/09_comparison/${marker}.vs_published.${rank}.name_translation.tsv", checkIfExists: true)) }
            COMPARE_ARMS(ASSIGN_ASV.out.taxa.combine(pub, by: 0), REFERENCE.out.db)
        }

        // Phase 6: the same analysis on arm 2's consensus table, on arm 1's, and on the published one.
        // the samplesheet, whichever way it was obtained: it carries seq_run, and without it the
        // sequencing batch is invisible to the ecology and the confound check cannot fire
        ECOLOGY_METADATA(file(params.published, checkIfExists: true), sheet)
        PUBLISHED_TABLE(file(params.published, checkIfExists: true))
        tables = CONSENSUS.out.table
                     .mix(PUBLISHED_TABLE.out.map { t -> tuple('published', t) })
        ECOLOGY(tables, ECOLOGY_METADATA.out)
    }
}
