#!/usr/bin/env nextflow

include { REFORMAT_MSA as REFORMAT_SEED_MSA     } from "./reformat_msa.nf"
include { REFORMAT_MSA as REFORMAT_HMMALIGN_MSA } from "./reformat_msa.nf"
include { ANNOTATE_MODELS                       } from "./annotate_models.nf"
include { PREDICT_STRUCTURES                    } from "./predict_structures.nf"
include { ANNOTATE_STRUCTURES                   } from "./annotate_structures.nf"

workflow ANNOTATE_FAMILIES {
    take:
    seed_msa_sto
    msa_sto

    main:
    seed_msa_sto
        .map { files ->
            String filePath = files[0]
            int lastIndex = filePath.lastIndexOf('/')
            String seed_msa_dir = filePath.substring(0, lastIndex + 1)
            return [ [id:"seed_msa"], file(seed_msa_dir) ]
        }
        .set { seed_msa_ch }
    
    msa_sto
        .map { files ->
            String filePath = files[0]
            int lastIndex = filePath.lastIndexOf('/')
            String msa_dir = filePath.substring(0, lastIndex + 1)
            return [ [id:"msa"], file(msa_dir) ]
        }
        .set { hmmalign_msa_ch }

    fa_seed_msa_ch = REFORMAT_SEED_MSA(seed_msa_ch).fa_ch
    fa_msa_ch      = REFORMAT_HMMALIGN_MSA(hmmalign_msa_ch).fa_ch
    pfam_hits      = ANNOTATE_MODELS(fa_seed_msa_ch)
    structure_ch   = PREDICT_STRUCTURES(hmmalign_msa_ch)
    pdb_ch         = structure_ch.pdb_ch
    scores_ch      = structure_ch.scores_ch
    cif_ch         = structure_ch.cif_ch
    foldseek_hits  = ANNOTATE_STRUCTURES(pdb_ch)

    emit:
    pfam_hits
    foldseek_hits
    scores_ch
    cif_ch
    fa_seed_msa_ch
    fa_msa_ch
}
