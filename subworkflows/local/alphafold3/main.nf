include { PREPARE_ALPHAFOLD3_DBS   } from '../../../subworkflows/local/prepare_alphafold3_dbs'
include { FASTA_TO_ALPHAFOLD3_JSON } from '../../../modules/local/fasta_to_alphafold3_json'
include { RUN_ALPHAFOLD3           } from '../../../modules/local/run_alphafold3'
include { MMCIF2PDB                } from '../../../modules/local/mmcif2pdb/main.nf'

workflow ALPHAFOLD3 {
    take:
    fasta
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

    PREPARE_ALPHAFOLD3_DBS (
        alphafold3_db,
        alphafold3_params_path,
        alphafold3_small_bfd_path,
        alphafold3_mgnify_path,
        alphafold3_pdb_mmcif_path,
        alphafold3_uniref90_path,
        alphafold3_pdb_seqres_path,
        alphafold3_uniprot_path,
        alphafold3_small_bfd_link,
        alphafold3_mgnify_link,
        alphafold3_pdb_mmcif_link,
        alphafold3_uniref90_link,
        alphafold3_pdb_seqres_link,
        uniprot_link
    )
    ch_versions = ch_versions.mix( PREPARE_ALPHAFOLD3_DBS.out.versions )

    FASTA_TO_ALPHAFOLD3_JSON( fasta )
    ch_versions = ch_versions.mix( FASTA_TO_ALPHAFOLD3_JSON.out.versions )

    RUN_ALPHAFOLD3( FASTA_TO_ALPHAFOLD3_JSON.out.json,
        PREPARE_ALPHAFOLD3_DBS.out.params,
        PREPARE_ALPHAFOLD3_DBS.out.small_bfd,
        PREPARE_ALPHAFOLD3_DBS.out.mgnify,
        PREPARE_ALPHAFOLD3_DBS.out.pdb_mmcif,
        PREPARE_ALPHAFOLD3_DBS.out.uniref90,
        PREPARE_ALPHAFOLD3_DBS.out.pdb_seqres,
        PREPARE_ALPHAFOLD3_DBS.out.uniprot )
    ch_versions = ch_versions.mix( RUN_ALPHAFOLD3.out.versions )

    MMCIF2PDB( RUN_ALPHAFOLD3.out.top_ranked_cif )
    ch_versions = ch_versions.mix( MMCIF2PDB.out.versions )

    emit:
    versions = ch_versions
    // scores   = ch_scores
    pdb      = MMCIF2PDB.out.pdb
    cif      = RUN_ALPHAFOLD3.out.top_ranked_cif
}
