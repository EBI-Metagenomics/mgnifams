include { ESMFOLD    } from '../../../subworkflows/local/esmfold'
// include { ALPHAFOLD3 } from '../../../subworkflows/local/alphafold3'

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

    ESMFOLD( ch_fasta, pdb_chunk_size, esmfold_db, esmfold_params_path, esmfold_3B_v1, \
        esm2_t36_3B_UR50D, esm2_t36_3B_UR50D_contact_regression, num_recycles_esmfold, \
        pdb_chunk_size_long, outdir)
    ch_versions = ch_versions.mix( ESMFOLD.out.versions )
    ch_scores = ESMFOLD.out.scores
    ch_pdb = ESMFOLD.out.pdb
    ch_cif = ESMFOLD.out.cif

    emit:
    versions = ch_versions
    scores   = ch_scores
    pdb      = ch_pdb
    cif      = ch_cif
}
