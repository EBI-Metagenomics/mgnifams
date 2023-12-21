#!/usr/bin/env nextflow

include { EXTRACT_FIRST_STOCKHOLM_SEQUENCES } from "$launchDir/modules/general.nf"
include { ESMFOLD                           } from "$launchDir/modules/esmfold/main.nf"
include { PARSE_CIF                         } from "$launchDir/modules/esmfold/parse_cif.nf"
include { FOLDCOMP_COMPRESS                 } from "$launchDir/modules/foldcomp/compress/main.nf"

workflow PREDICT_STRUCTURES {
    take:
    msa
    unannotated_ids
    mode
    
    main:
    family_reps_fasta = EXTRACT_FIRST_STOCKHOLM_SEQUENCES(msa, unannotated_ids, mode)

    pdb_fasta_chunks_ch = family_reps_fasta.splitFasta( by: params.pdb_chunk_size, file: true )
    pdb_fasta_chunks_ch
        .map { filepath ->
            def parts = filepath.baseName.split('\\.')
            def number = parts[1]
            def id = "pdb${number}"
            return [ [id:id], [file(filepath)] ]
        }
        .set { input }
    pdb_ch = ESMFOLD(input).pdb
    cif_ch = PARSE_CIF(pdb_ch).cif
    fcz_ch = FOLDCOMP_COMPRESS(pdb_ch).fcz

    emit:
    pdb_ch
    cif_ch
    fcz_ch
}
