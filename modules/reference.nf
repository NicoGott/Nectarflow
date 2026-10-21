/*
 * The reference library, shared by both arms so a comparison between them is a comparison of
 * methods and not of libraries. Process text is unchanged from main.nf, so the work cache is shared.
 */
nextflow.enable.dsl = 2

process REFERENCE {
    // Built from the pinned manifest: Barcode UK for the UK flora, GenBank for the rest of the
    // Plant Atlas 2020. Nothing here is searched or fetched,
    // so the library is the same every time and the same for both arms. See docs/replica.md.
    publishDir "${params.outdir}/06_reference", mode: 'copy'
    conda "conda-forge::python=3.11 ${params.conda_blast}"
    input:  path manifest
            path sequences
    output: path 'refdb', emit: db
    script:
    def markers = params.markers.keySet().join(',')
    """
    build_reference.py --manifest ${manifest} --vendored ${sequences}
    mkdir refdb
    for m in \$(echo ${markers} | tr ',' ' '); do
      makeblastdb -in \$m.fasta -dbtype nucl -out refdb/\$m
      mv \$m.fasta \$m.taxonomy.tsv refdb/
    done
    cp ${manifest} refdb/reference_manifest.tsv
    """
    stub: "mkdir refdb; for m in ${params.markers.keySet().join(' ')}; do touch refdb/\$m.taxonomy.tsv; done"
}
