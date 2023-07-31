#!/usr/bin/env nextflow

include { HMMSCAN } from "$baseDir/modules/hmm.nf"
include { INTERPROSCAN } from "$baseDir/modules/interproscan.nf"
include { EGGNOG_MAPPER } from "$baseDir/modules/eggnog.nf"
// TODO protENN
include { EXPORT_INTERPRO_ANNOTATIONS_CSV; EXPORT_EGGNOG_ANNOTATIONS_CSV; EXPORT_UNIPROT_ANNOTATIONS_CSV; CONCAT_ANNOTATIONS } from "$baseDir/modules/exporting.nf"

process getTopDebugLines { // TODO remove in Codon
    input:
    path(fasta)

    output:
    path("top_reps.fasta")

    script:
    lines = params.debug_top_n_reps * 2
    """
    head -n ${lines} ${fasta} > top_reps.fasta
    """
}

workflow annotate_unknown {
    take:
    unknown_reps_fasta
    build_ch
    
    main:
    Channel
        .fromPath(params.uniprot_sprot_fasta_path) 
        .set { uniprot_sprot_fasta_path }

    top_reps = getTopDebugLines(unknown_reps_fasta) // TODO remove
    interpro_ch = INTERPROSCAN(top_reps) // TODO: reps_fa
    eggnog_ch = EGGNOG_MAPPER(top_reps) // TODO:  reps_fa
    tblout_ch = HMMSCAN(build_ch, uniprot_sprot_fasta_path.first()).tblout_ch
    interpro_csv = EXPORT_INTERPRO_ANNOTATIONS_CSV(interpro_ch)
    eggnog_csv = EXPORT_EGGNOG_ANNOTATIONS_CSV(eggnog_ch)
    uniprot_csv = EXPORT_UNIPROT_ANNOTATIONS_CSV(tblout_ch)
    CONCAT_ANNOTATIONS(interpro_csv.concat(eggnog_csv, uniprot_csv.collect()).collect(), 'unknown')
}