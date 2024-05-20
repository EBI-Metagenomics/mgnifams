#!/usr/bin/env nextflow

include { EXTRACT_FIRST_STOCKHOLM_SEQUENCES } from "${params.moduleDir}/family/main.nf"
include { ESMFOLD                           } from "${params.moduleDir}/esmfold/main.nf"
include { EXTRACT_ESMFOLD_SCORES            } from "${params.moduleDir}/esmfold/extract_esmfold_scores.nf"
include { PARSE_CIF                         } from "${params.moduleDir}/esmfold/parse_cif.nf"

workflow PREDICT_STRUCTURES {
    take:
    msa_sto_ch
    
    main:
    family_reps_fa_ch = EXTRACT_FIRST_STOCKHOLM_SEQUENCES(msa_sto_ch)
    fasta_chunks_ch = family_reps_fa_ch.splitFasta( by: params.pdb_chunk_size, file: true )
    fasta_chunks_ch
        .map { filepath ->
            def parts = filepath.baseName.split('\\.')
            def number = parts[1]
            def id = "pdb${number}"
            return [ [id:id], [file(filepath)] ]
        }
        .set { pdb_ch }
    esmfold_result = ESMFOLD(pdb_ch, params.compute_mode)
    pdb_ch = esmfold_result.pdb
    scores_ch = EXTRACT_ESMFOLD_SCORES(esmfold_result.scores).csv.map { meta, filepath -> filepath }
    scores_ch.collectFile(name: "pdb_scores.csv", storeDir: params.outDir + "/structures")
    PARSE_CIF(pdb_ch)

    emit:
    pdb_ch
}
