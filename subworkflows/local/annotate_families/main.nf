#!/usr/bin/env nextflow

include { REFORMAT_MSA as REFORMAT_SEED_MSA     } from "../../../subworkflows/local/reformat_msa"
include { REFORMAT_MSA as REFORMAT_HMMALIGN_MSA } from "../../../subworkflows/local/reformat_msa"
include { ANNOTATE_MODELS                       } from "../../../subworkflows/local/annotate_models"
include { PREDICT_STRUCTURES                    } from "../../../subworkflows/local/predict_structures"
include { ANNOTATE_STRUCTURES                   } from "../../../subworkflows/local/annotate_structures"

workflow ANNOTATE_FAMILIES {
    take:
    seed_msa_sto
    msa_sto

    main:
    ch_versions = Channel.empty()

    seed_msa_sto
        .map { meta, files ->
            String filePath = files[0]
            int lastIndex = filePath.lastIndexOf('/')
            String seed_msa_dir = filePath.substring(0, lastIndex + 1)
            [ [id:"seed_msa"], file(seed_msa_dir) ]
        }
        .set { seed_msa_ch }
    
    msa_sto
        .map { meta, files ->
            String filePath = files[0]
            int lastIndex = filePath.lastIndexOf('/')
            String msa_dir = filePath.substring(0, lastIndex + 1)
            [ [id:"msa"], file(msa_dir) ]
        }
        .set { hmmalign_msa_ch }

    ch_fa_seed_msa = REFORMAT_SEED_MSA(seed_msa_ch).fasta
    ch_versions = ch_versions.mix( REFORMAT_SEED_MSA.out.versions )

    ch_fa_msa = REFORMAT_HMMALIGN_MSA(hmmalign_msa_ch).fasta
    ch_versions = ch_versions.mix( REFORMAT_HMMALIGN_MSA.out.versions )

    ch_pfam_hits = ANNOTATE_MODELS(ch_fa_seed_msa).pfam_hits
    ch_versions = ch_versions.mix( ANNOTATE_MODELS.out.versions )
    
    ch_structures = PREDICT_STRUCTURES(hmmalign_msa_ch)
    ch_versions = ch_versions.mix( PREDICT_STRUCTURES.out.versions )

    ch_pdb    = ch_structures.pdb
    ch_scores = ch_structures.scores
    ch_cif    = ch_structures.cif

    ch_foldseek_hits = ANNOTATE_STRUCTURES(ch_pdb).foldseek_hits
    ch_versions = ch_versions.mix( ANNOTATE_STRUCTURES.out.versions )

    emit:
    versions      = ch_versions
    pfam_hits     = ch_pfam_hits
    foldseek_hits = ch_foldseek_hits
    scores        = ch_scores
    cif           = ch_cif
    fa_seed_msa   = ch_fa_seed_msa
    fa_msa        = ch_fa_msa
}
