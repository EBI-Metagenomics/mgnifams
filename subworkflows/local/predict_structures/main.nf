include { ESMFOLD    } from '../../../subworkflows/local/esmfold'
include { ALPHAFOLD3 } from '../../../subworkflows/local/alphafold3'

workflow PREDICT_STRUCTURES {
    take:
    fasta
    fold_mode
    pdb_chunk_size
    esmfold_db
    esmfold_params_path
    esmfold_3B_v1
    esm2_t36_3B_UR50D
    esm2_t36_3B_UR50D_contact_regression
    num_recycles_esmfold
    pdb_chunk_size_long
    outdir
    alphafold3_db
    alphafold3_params_path
    alphafold3_small_bfd_path
    alphafold3_mgnify_path
    alphafold3_pdb_mmcif_path
    alphafold3_uniref90_path
    alphafold3_pdb_seqres_path
    alphafold3_uniprot_path
    alphafold3_small_bfd_link
    alphafold3_mgnify_link
    alphafold3_pdb_mmcif_link
    alphafold3_uniref90_link
    alphafold3_pdb_seqres_link
    uniprot_link
    
    main:
    ch_versions = Channel.empty()

    if (fold_mode == 'alphafold3') {
        pdb_chunk_size = 1 // can only predict one protein at a time with alphafold3
    }

    ch_fasta = fasta
        .map { meta, file_path ->
            [ meta, file_path.splitFasta( by: pdb_chunk_size, file: true ) ]
        }
        .transpose()
        .map { meta, file_path ->
            [ [id: meta.id, chunk: file(file_path, checkIfExists: true).getBaseName().split('\\.')[-1]], file_path ]
        }

    if (fold_mode == 'esmfold') {
        ESMFOLD( ch_fasta, pdb_chunk_size, esmfold_db, esmfold_params_path, esmfold_3B_v1, \
        esm2_t36_3B_UR50D, esm2_t36_3B_UR50D_contact_regression, num_recycles_esmfold, \
        pdb_chunk_size_long, outdir)
        ch_versions = ch_versions.mix( ESMFOLD.out.versions )
        ch_scores = ESMFOLD.out.scores
        ch_pdb = ESMFOLD.out.pdb
        ch_cif = ESMFOLD.out.cif
    } else if (fold_mode == 'alphafold3') {
        ALPHAFOLD3( ch_fasta, alphafold3_db, alphafold3_params_path, alphafold3_small_bfd_path, \
        alphafold3_mgnify_path, alphafold3_pdb_mmcif_path, alphafold3_uniref90_path, alphafold3_pdb_seqres_path, \
        alphafold3_uniprot_path, alphafold3_small_bfd_link, alphafold3_mgnify_link, alphafold3_pdb_mmcif_link, \
        alphafold3_uniref90_link, alphafold3_pdb_seqres_link, uniprot_link)
        ch_versions = ch_versions.mix( ALPHAFOLD3.out.versions )
        ch_pdb = ALPHAFOLD3.out.pdb
        ch_cif = ALPHAFOLD3.out.cif
    }

    emit:
    versions = ch_versions
    // scores   = ch_scores
    pdb      = ch_pdb
    cif      = ch_cif
}
