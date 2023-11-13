#!/usr/bin/env nextflow

include { FIND_FASTA_BY_ID } from "$baseDir/modules/general.nf"
include { ESMFOLD } from "$baseDir/modules/esmfold/main.nf"

workflow annotate_structures {
    take:
    ids
    fasta
    
    main:
    FIND_FASTA_BY_ID(ids, fasta)
    pdb_fasta_chunks_ch = FIND_FASTA_BY_ID.out.splitFasta( by: params.pdb_chunk_size, file: true )
    pdb_fasta_chunks_ch
        .map { filepath ->
            def parts = filepath.baseName.split('\\.')
            def number = parts[1]
            def id = "pdb${number}"
            return [ [id:id], [file(filepath)] ]
        }
        .set { input }
    ESMFOLD(input)
    // esm_output = ESMFOLD.map { meta, pdb -> pdb }.view()
    // esm_output.view()

    // emit:
    // esm_output
}