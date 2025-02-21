#!/usr/bin/env nextflow

include { MMSEQS_CREATEDB  } from "../../../modules/nf-core/mmseqs/createdb/main"
include { MMSEQS_LINCLUST  } from "../../../modules/nf-core/mmseqs/linclust/main"
include { MMSEQS_CREATETSV } from "../../../modules/nf-core/mmseqs/createtsv/main"

workflow EXECUTE_CLUSTERING {
    take:
    mgnifams_input_fasta

    main:
    ch_versions = Channel.empty()

    MMSEQS_CREATEDB( mgnifams_input_fasta )
    ch_versions = ch_versions.mix( MMSEQS_CREATEDB.out.versions )

    MMSEQS_LINCLUST( MMSEQS_CREATEDB.out.db )
    ch_versions = ch_versions.mix( MMSEQS_LINCLUST.out.versions )

    MMSEQS_CREATETSV( MMSEQS_LINCLUST.out.db_cluster, MMSEQS_CREATEDB.out.db, MMSEQS_CREATEDB.out.db )
    ch_versions = ch_versions.mix( MMSEQS_CREATETSV.out.versions )
    
    emit:
    versions      = ch_versions
    clusters_tsv  = MMSEQS_CREATETSV.out.tsv
}
