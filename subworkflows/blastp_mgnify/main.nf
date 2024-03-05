#!/usr/bin/env nextflow

include { DIAMOND_MAKEDB } from "$launchDir/modules/diamond/makedb/main"
include { DIAMOND_BLASTP } from "$launchDir/modules/diamond/blastp/main"
include { MATCH_CLUSTERS } from "$launchDir/modules/match_clusters.nf"

workflow BLASTP_MGNIFY {
    take:
    anti_defence_proteins
    mgnifams_input
    clusters_tsv

    main:
    DIAMOND_MAKEDB ( [ [ id:'mgnifams' ], mgnifams_input ], [], [], [] )
    blast_columns = '' // defaults: qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore
    DIAMOND_BLASTP ( [ [ id:'anti_defence_proteins' ], anti_defence_proteins ], DIAMOND_MAKEDB.out.db, 'txt', blast_columns )
    blastp_tsv = DIAMOND_BLASTP.out.txt

    matched_clusters = MATCH_CLUSTERS(blastp_tsv, clusters_tsv).csv

    emit:
    matched_clusters
}
