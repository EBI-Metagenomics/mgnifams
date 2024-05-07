#!/usr/bin/env nextflow

include { EXTRACT_FIRST_STOCKHOLM_SEQUENCES } from "$launchDir/modules/family/main.nf"
include { ESMFOLD                           } from "$launchDir/modules/esmfold/main.nf"
include { EXTRACT_ESMFOLD_SCORES            } from "$launchDir/modules/esmfold/extract_esmfold_scores.nf"
include { PARSE_CIF                         } from "$launchDir/modules/esmfold/parse_cif.nf"

workflow PREDICT_STRUCTURES {
    take:
    msa_ch
    
    main:
    family_reps_fasta = EXTRACT_FIRST_STOCKHOLM_SEQUENCES(msa_ch)

    pdb_fasta_chunks_ch = family_reps_fasta.splitFasta( by: params.pdb_chunk_size, file: true )
    pdb_fasta_chunks_ch
        .map { filepath ->
            def parts = filepath.baseName.split('\\.')
            def number = parts[1]
            def id = "pdb${number}"
            return [ [id:id], [file(filepath)] ]
        }
        .set { input }
    esmfold_result = ESMFOLD(input, params.compute_mode)
    pdb_ch = esmfold_result.pdb
    EXTRACT_ESMFOLD_SCORES(esmfold_result.scores)
    cif_ch = PARSE_CIF(pdb_ch).cif

    emit:
    pdb_ch
    cif_ch
}
