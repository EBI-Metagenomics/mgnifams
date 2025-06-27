include { PREPARE_ESMFOLD_DBS } from '../../../subworkflows/local/prepare_esmfold_dbs'
// include { ESMFOLD             } from '../../../subworkflows/local/esmfold'

// include { RUN_ESMFOLD_CONDA } from '../../../modules/local/run_esmfold_conda/main'
// include { RUN_ESMFOLD       } from '../../../modules/local/run_esmfold/main'

// include { EXTRACT_LONG_FA        } from '../../../modules/local/extract_long_fa/main'
// include { ESMFOLD_CPU            } from '../../../modules/local/esmfold/esmfold_cpu/main'
// include { EXTRACT_ESMFOLD_SCORES } from '../../../modules/local/extract_esmfold_scores/main'
// include { PARSE_CIF              } from '../../../modules/local/parse_cif/main'

workflow PREDICT_STRUCTURES {
    take:
    fasta
    pdb_chunk_size

    esmfold_db
    esmfold_params_path
    esmfold_3B_v1
    esm2_t36_3B_UR50D
    esm2_t36_3B_UR50D_contact_regression

    // compute_mode
    // pdb_chunk_size_long // TODO maybe remove
    // outdir
    
    main:
    ch_versions = Channel.empty()

    ch_fasta = fasta
        .map { meta, file_path ->
            [ meta, file_path.splitFasta( by: pdb_chunk_size, file: true ) ]
        }
        .transpose()
        .map { meta, file_path ->
            [ [id: meta.id, chunk: file(file_path, checkIfExists: true).getBaseName().split('\\.')[-1]], file_path ]
        }

    PREPARE_ESMFOLD_DBS( esmfold_db, esmfold_params_path, esmfold_3B_v1, \
        esm2_t36_3B_UR50D, esm2_t36_3B_UR50D_contact_regression )
    ch_versions = ch_versions.mix( PREPARE_ESMFOLD_DBS.out.versions )

    // ESMFOLD (
    //     ch_fasta,
    //     ch_versions,
    //     PREPARE_ESMFOLD_DBS.out.params,
    //     num_recycles_esmfold
    // )
    // ch_versions = ch_versions.mix( ESMFOLD.out.versions )

    // if (workflow.profile.contains("conda")) {
    //     esmfold_result = RUN_ESMFOLD_CONDA( ch_fasta, compute_mode )
    //     ch_versions = ch_versions.mix( RUN_ESMFOLD_CONDA.out.versions )
    // } else {
    //     esmfold_result = RUN_ESMFOLD( ch_fasta, compute_mode )
    //     ch_versions = ch_versions.mix( RUN_ESMFOLD.out.versions )
    // }
    
    // // Long sequences that cannot be run on GPU
    // fa_ch
    //     .map { meta, files ->
    //         files
    //     }
    //     .collect()
    //     .map { file ->
    //         [ [id:"esm_scores"], file ]
    //     }
    //     .set { fa_paths_ch }

    // esmfold_result.scores
    //     .map { meta, files ->
    //         files
    //     }
    //     .collect()
    //     .map { file ->
    //         [ [id:"esm_scores"], file ]
    //     }
    //     .set { score_paths_ch }
    
    // long_reps_fa_ch = EXTRACT_LONG_FA(fa_paths_ch, score_paths_ch)
    // // TODO ch_versions = ch_versions.mix( EXTRACT_LONG_FA.out.versions )

    // long_reps_fa_ch
    //     .map { meta, filepath ->
    //         [ meta, filepath.splitFasta( by: pdb_chunk_size_long, file: true ) ]
    //     }
    //     .transpose()
    //     .map { meta, filepath ->
    //         [ [id: meta.id, chunk: filepath.getBaseName().split('\\.')[-1]], filepath ]
    //     }
    //     .set { fa_long_ch }

    // esmfold_long_result = ESMFOLD_CPU(fa_long_ch)
    // ch_versions = ch_versions.mix( ESMFOLD_CPU.out.versions )
    // // End long sequences

    // ch_scores = EXTRACT_ESMFOLD_SCORES(esmfold_result.scores.concat(esmfold_long_result.scores)).csv
    // // TODO ch_versions = ch_versions.mix( EXTRACT_ESMFOLD_SCORES.out.versions )
    
    // ch_scores = ch_scores
    //     .map { meta, file ->
    //         file
    //     }
    //     .collectFile(name: "pdb_scores.csv", storeDir: outdir + "/structures")
    //     .map { file ->
    //         [ [id: "scores"], file ]
    //     }
    
    // ch_pdb = esmfold_result.pdb.concat(esmfold_long_result.pdb)
    // ch_cif = PARSE_CIF(ch_pdb)
    // // TODO ch_versions = ch_versions.mix( PARSE_CIF.out.versions )

    // ch_cif = ch_cif
    //     .map { meta, filepath ->
    //         filepath }
    //     .collect()
    //     .map { file ->
    //         [ [id: "cif"], file ]
    //     }

    emit:
    versions = ch_versions
    // esm_multiqc = ESMFOLD.out.multiqc_report.collect()
    // pdb         = ESMFOLD.out.pdb

    // scores   = ch_scores
    // pdb      = ch_pdb
    // cif      = ch_cif
}
