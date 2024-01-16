#!/usr/bin/env nextflow

include { MMSEQS_CREATEDB  } from "$launchDir/modules/mmseqs/createdb/main.nf"
include { MMSEQS_LINCLUST  } from "$launchDir/modules/mmseqs/linclust/main.nf"
include { CREATETSV        } from "$launchDir/modules/mmseqs2.nf"

workflow EXECUTE_CLUSTERING {
    take:
    fasta_file

    main:
    mmseqs_db = MMSEQS_CREATEDB( fasta_file.map { fasta -> [ [ id:'mmseqs_db' ], fasta ] } ).db
    mmseqs_families = MMSEQS_LINCLUST(mmseqs_db).db_cluster
    families_tsv = CREATETSV(mmseqs_db, mmseqs_families)
    
    emit:
    families_tsv
}
