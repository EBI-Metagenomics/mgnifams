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
    ch_versions = Channel.empty()

    family_reps_fa_ch = EXTRACT_FIRST_STOCKHOLM_SEQUENCES_FROM_FOLDER(msa_sto_ch)
    // TODO ch_versions = ch_versions.mix( EXTRACT_FIRST_STOCKHOLM_SEQUENCES_FROM_FOLDER.out.versions )

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
    ch_versions = ch_versions.mix( ESMFOLD.out.versions )
    
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
    // TODO ch_versions = ch_versions.mix( EXTRACT_LONG_FA.out.versions )

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
    ch_versions = ch_versions.mix( ESMFOLD_CPU.out.versions )
    // End long sequences

    ch_scores = EXTRACT_ESMFOLD_SCORES(esmfold_result.scores.concat(esmfold_long_result.scores)).csv
    // TODO ch_versions = ch_versions.mix( EXTRACT_ESMFOLD_SCORES.out.versions )
    
    ch_scores = ch_scores
        .map { meta, file ->
            file
        }
        .collectFile(name: "pdb_scores.csv", storeDir: params.outdir + "/structures")
        .map { file ->
            [ [id: "scores"], file ]
        }
    
    ch_pdb = esmfold_result.pdb.concat(esmfold_long_result.pdb)
    ch_cif = PARSE_CIF(ch_pdb)
    // TODO ch_versions = ch_versions.mix( PARSE_CIF.out.versions )

    ch_cif = ch_cif
        .map { meta, filepath ->
            filepath }
        .collect()
        .map { file ->
            [ [id: "cif"], file ]
        }

    emit:
    versions = ch_versions
    scores   = ch_scores
    pdb      = ch_pdb
    cif      = ch_cif
}
