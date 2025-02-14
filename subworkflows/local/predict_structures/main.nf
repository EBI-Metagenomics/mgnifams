#!/usr/bin/env nextflow

include { EXTRACT_FIRST_STOCKHOLM_SEQUENCES_FROM_FOLDER } from "../../../modules/local/extract_first_stockholm_sequences_from_folder/main"
include { ESMFOLD                                       } from "../../../modules/local/esmfold/esmfold/main"
include { EXTRACT_LONG_FA                               } from "../../../modules/local/extract_long_fa/main"
include { ESMFOLD_CPU                                   } from "../../../modules/local/esmfold/esmfold_cpu/main"
include { EXTRACT_ESMFOLD_SCORES                        } from "../../../modules/local/extract_esmfold_scores/main"
include { PARSE_CIF                                     } from "../../../modules/local/parse_cif/main"

workflow PREDICT_STRUCTURES {
    take:
    msa_sto_ch
    
    main:
    family_reps_fa_ch = EXTRACT_FIRST_STOCKHOLM_SEQUENCES_FROM_FOLDER(msa_sto_ch)
    family_reps_fa_ch
        .map { meta, filepath ->
            [ meta, filepath.splitFasta( by: params.pdb_chunk_size, file: true ) ]
        }
        .transpose()
        .map { meta, filepath ->
            [ [id: meta.id, chunk: filepath.getBaseName().split('\\.')[-1]], filepath ]
        }
        .set { fa_ch }

    esmfold_result = ESMFOLD(fa_ch, params.compute_mode)
    
    // Long sequences that cannot be run on GPU
    fa_ch
        .map { meta, files ->
            files
        }
        .collect()
        .map { file ->
            [ [id:"esm_scores"], file ]
        }
        .set { fa_paths_ch }

    esmfold_result.scores
        .map { meta, files ->
            files
        }
        .collect()
        .map { file ->
            [ [id:"esm_scores"], file ]
        }
        .set { score_paths_ch }
    
    long_reps_fa_ch = EXTRACT_LONG_FA(fa_paths_ch, score_paths_ch)

    long_reps_fa_ch
        .map { meta, filepath ->
            [ meta, filepath.splitFasta( by: params.pdb_chunk_size_long, file: true ) ]
        }
        .transpose()
        .map { meta, filepath ->
            [ [id: meta.id, chunk: filepath.getBaseName().split('\\.')[-1]], filepath ]
        }
        .set { fa_long_ch }

    esmfold_long_result = ESMFOLD_CPU(fa_long_ch)
    // End long sequences

    scores_ch = EXTRACT_ESMFOLD_SCORES(esmfold_result.scores.concat(esmfold_long_result.scores)).csv
    
    scores_ch
        .map { meta, file ->
            file
        }
        .collectFile(name: "pdb_scores.csv", storeDir: params.outdir + "/structures")
        .map { file ->
            [ [id: "scores"], file ]
        }
        .set { scores_ch }
    
    pdb_ch = esmfold_result.pdb.concat(esmfold_long_result.pdb)
    cif_ch = PARSE_CIF(pdb_ch)
    cif_ch
        .map { meta, filepath ->
            filepath }
        .collect()
        .map { file ->
            [ [id: "cif"], file ]
        }
        .set { cif_ch }

    emit:
    scores_ch
    pdb_ch
    cif_ch
}
