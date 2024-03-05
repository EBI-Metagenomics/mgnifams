#!/usr/bin/env nextflow

include { DIAMOND_MAKEDB } from "$launchDir/modules/diamond/makedb/main"
include { DIAMOND_BLASTP } from "$launchDir/modules/diamond/blastp/main"

workflow BLASTP_MGNIFY {
    take:
    anti_defence_proteins
    mgnifams_input

    main:
    DIAMOND_MAKEDB ( [ [ id:'mgnifams' ], mgnifams_input ], [], [], [] )
    blast_columns = '' // defaults: qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore
    DIAMOND_BLASTP ( [ [ id:'anti_defence_proteins' ], anti_defence_proteins ], DIAMOND_MAKEDB.out.db, 'txt', blast_columns )
    blastp_tsv = DIAMOND_BLASTP.out.txt

    emit:
    blastp_tsv
}
