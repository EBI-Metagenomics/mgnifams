include { PREPARE_ESMFOLD_DBS            } from '../../../subworkflows/local/prepare_esmfold_dbs'
include { RUN_ESMFOLD                    } from '../../../modules/local/run_esmfold'
include { EXTRACT_CUDA_FAILED            } from '../../../modules/local/extract_cuda_failed/main'
include { RUN_ESMFOLD as RUN_ESMFOLD_CPU } from '../../../modules/local/run_esmfold'
include { EXTRACT_ESMFOLD_SCORES         } from '../../../modules/local/extract_esmfold_scores/main'
include { PARSE_CIF                      } from '../../../modules/local/parse_cif/main'

workflow PREDICT_STRUCTURES {
    take:
    fasta
    pdb_chunk_size
    esmfold_db
    esmfold_params_path
    esmfold_3B_v1
    esm2_t36_3B_UR50D
    esm2_t36_3B_UR50D_contact_regression
    num_recycles_esmfold
    pdb_chunk_size_long
    outdir
    
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

    RUN_ESMFOLD( ch_fasta, PREPARE_ESMFOLD_DBS.out.params, num_recycles_esmfold )
    ch_versions = ch_versions.mix( RUN_ESMFOLD.out.versions )

    // Identify CUDA failed very long sequences, and run on CPU
    ch_scores = RUN_ESMFOLD.out.scores
        .map { meta, files ->
            files
        }
        .collect()
        .map { file ->
            [ [id:"esm_scores"], file ]
        }
    
    EXTRACT_CUDA_FAILED(fasta, ch_scores)
    ch_versions = ch_versions.mix( EXTRACT_CUDA_FAILED.out.versions )

    ch_fasta_long = EXTRACT_CUDA_FAILED.out.fasta
        .map { meta, file_path ->
            [ meta, file_path.splitFasta( by: pdb_chunk_size_long, file: true ) ]
        }
        .transpose()
        .map { meta, file_path ->
            [ [id: meta.id, chunk: file(file_path, checkIfExists: true).getBaseName().split('\\.')[-1]], file_path ]
        }

    RUN_ESMFOLD_CPU( ch_fasta_long, PREPARE_ESMFOLD_DBS.out.params, num_recycles_esmfold )
    ch_versions = ch_versions.mix( RUN_ESMFOLD_CPU.out.versions )

    EXTRACT_ESMFOLD_SCORES( RUN_ESMFOLD.out.scores.concat(RUN_ESMFOLD_CPU.out.scores) )
    ch_versions = ch_versions.mix( EXTRACT_ESMFOLD_SCORES.out.versions )
    
    ch_scores = EXTRACT_ESMFOLD_SCORES.out.csv
        .map { meta, file ->
            file
        }
        .collectFile(name: "pdb_scores.csv", storeDir: outdir + "/structures/esmfold/", keepHeader: true)
        .map { file ->
            [ [id: "scores"], file ]
        }

    PARSE_CIF( RUN_ESMFOLD.out.pdb.concat(RUN_ESMFOLD_CPU.out.pdb) )
    ch_versions = ch_versions.mix( PARSE_CIF.out.versions )

    ch_pdb = RUN_ESMFOLD.out.pdb.concat(RUN_ESMFOLD_CPU.out.pdb)
        .map { meta, file_path ->
            file_path }
        .collect()
        .map { file ->
            [ [id: "pdb"], file ]
        }

    ch_cif = PARSE_CIF.out.cif
        .map { meta, file_path ->
            file_path }
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
