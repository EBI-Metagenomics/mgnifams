#!/usr/bin/env nextflow

include { MMSEQS_CREATEDB  } from "../modules/nf-core/mmseqs/createdb/main.nf"
include { MMSEQS_LINCLUST  } from "../modules/nf-core/mmseqs/linclust/main.nf"
include { MMSEQS_CREATETSV } from "../modules/nf-core/mmseqs/createtsv/main.nf"

workflow EXECUTE_CLUSTERING {
    take:
    mgnifams_input_fasta

    main:
    MMSEQS_CREATEDB( mgnifams_input_fasta )
    MMSEQS_LINCLUST( MMSEQS_CREATEDB.out.db )
    MMSEQS_CREATETSV( MMSEQS_LINCLUST.out.db_cluster, MMSEQS_CREATEDB.out.db, MMSEQS_CREATEDB.out.db )
    
    emit:
    clusters_tsv  = MMSEQS_CREATETSV.out.tsv
    num_sequences = MMSEQS_CREATETSV.out.num_sequences
}
