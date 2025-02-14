#!/usr/bin/env nextflow

include { GENERATE_FAMILIES_PARALLEL } from "../../../subworkflows/local/generate_families_parallel"
include { FLAG_TRANSMEMBRANE         } from "../../../subworkflows/local/flag_transmembrane"
include { MOVE_TO_DIR                } from "../../../modules/local/move_to_dir/main"
include { REMOVE_REDUNDANCY          } from "../../../subworkflows/local/remove_redundancy"

workflow GENERATE_NONREDUNDANT_FAMILIES {
    take:
    clusters_tsv
    checked_clusters
    mgnifams_input_fa

    main:
    families_ch = GENERATE_FAMILIES_PARALLEL(clusters_tsv, [], mgnifams_input_fa)

    families_ch.msa_sto
        .map { meta, files ->
            files
        }
        .collect()
        .map { file ->
            [ [id:"msa_sto"], file ]
        }
        .set { msa_sto_ch }
    
    tm_ch     = FLAG_TRANSMEMBRANE(msa_sto_ch)
    rep_fa_ch = tm_ch.fa_ch
    tm_ids_ch = tm_ch.tm_ids_ch
    prob_ids  = tm_ch.prob_ids
    
    families_ch.seed_msa_sto
        .map { meta, files ->
            files
        }
        .collect()
        .map { file ->
            [ [id:"seed_msa_sto"], file ]
        }
        .set { seed_msa_sto_ch }

    seed_msa_sto_dir = MOVE_TO_DIR(seed_msa_sto_ch, "seed_msa_sto")

    families_ch.hmm
        .map { meta, files ->
            files
        }
        .collect()
        .map { file ->
            [ [id:"hmm"], file ]
        }
        .set { hmm_ch }

    families_ch.rf
        .map { meta, files ->
            files
        }
        .collect()
        .map { file ->
            [ [id:"rf"], file ]
        }
        .set { rf_ch }

    families_ch.domtblout
        .map { meta, files ->
            files
        }
        .collect()
        .map { file ->
            [ [id:"domtblout"], file ]
        }
        .set { domtblout_ch }

    families_ch.tsv
        .map { meta, files ->
            files
        }
        .collect()
        .map { file ->
            [ [id:"tsv"], file ]
        }
        .set { tsv_ch }

    families_ch.discarded
        .map { meta, files ->
            files
        }
        .collect()
        .map { file ->
            [ [id:"discarded"], file ]
        }
        .set { discarded_ch }

    families_ch.successful
        .map { meta, files ->
            files
        }
        .collect()
        .map { file ->
            [ [id:"successful"], file ]
        }
        .set { successful_ch }

    families_ch.converged
        .map { meta, files ->
            files
        }
        .collect()
        .map { file ->
            [ [id:"converged"], file ]
        }
        .set { converged_ch }

    families_ch.metadata
        .map { meta, files ->
            files
        }
        .collect()
        .map { file ->
            [ [id:"metadata"], file ]
        }
        .set { metadata_ch }

    families_ch.logs
        .map { meta, files ->
            files
        }
        .collect()
        .map { file ->
            [ [id:"logs"], file ]
        }
        .set { logs_ch }

    generated_families = REMOVE_REDUNDANCY(seed_msa_sto_dir, seed_msa_sto_ch, \
        msa_sto_ch, hmm_ch, \
        rf_ch, domtblout_ch, tsv_ch, \
        discarded_ch, successful_ch, converged_ch, \
        metadata_ch, logs_ch, tm_ids_ch, prob_ids, rep_fa_ch)

    emit:
    seed_msa_sto = generated_families.seed_msa_sto
    msa_sto      = generated_families.msa_sto
    metadata     = generated_families.metadata
    converged    = generated_families.converged
    tsv          = generated_families.tsv
    hmm          = generated_families.hmm
    rf           = generated_families.rf
}
