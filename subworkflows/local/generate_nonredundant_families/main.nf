#!/usr/bin/env nextflow

include { REFINE_FAMILIES_PARALLEL } from "../../../modules/local/refine_families_parallel/main" // TODO rename to GENERATE_FAMILIES
include { FLAG_TRANSMEMBRANE       } from "../../../subworkflows/local/flag_transmembrane"
include { MOVE_TO_DIR              } from "../../../modules/local/move_to_dir/main"
include { REMOVE_REDUNDANCY        } from "../../../subworkflows/local/remove_redundancy"

workflow GENERATE_NONREDUNDANT_FAMILIES {
    take:
    cluster_chunks
    mgnifams_fa

    main:
    ch_versions = Channel.empty()

    ch_families = REFINE_FAMILIES_PARALLEL( cluster_chunks, mgnifams_fa.first() )
    ch_versions = ch_versions.mix( REFINE_FAMILIES_PARALLEL.out.versions )

    // TODO multiMap?
    ch_msa_sto = ch_families.msa_sto
        .map { meta, files ->
            files
        }
        .collect()
        .map { file ->
            [ [id:"msa_sto"], file ]
        }
    
    ch_tm = FLAG_TRANSMEMBRANE( ch_msa_sto )
    ch_versions = ch_versions.mix( FLAG_TRANSMEMBRANE.out.versions )

    ch_rep_fa = ch_tm.fasta
    ch_tm_ids = ch_tm.tm_ids
    prob_ids = ch_tm.prob_ids
    
    ch_seed_msa_sto = ch_families.seed_msa_sto
        .map { meta, files ->
            files
        }
        .collect()
        .map { file ->
            [ [id:"seed_msa_sto"], file ]
        }

    seed_msa_sto_dir = MOVE_TO_DIR( ch_seed_msa_sto, "seed_msa_sto" )
    // TODO ch_versions = ch_versions.mix( MOVE_TO_DIR.out.versions )

    ch_hmm = ch_families.hmm
        .map { meta, files ->
            files
        }
        .collect()
        .map { file ->
            [ [id:"hmm"], file ]
        }

    ch_rf = ch_families.rf
        .map { meta, files ->
            files
        }
        .collect()
        .map { file ->
            [ [id:"rf"], file ]
        }

    ch_domtblout = ch_families.domtblout
        .map { meta, files ->
            files
        }
        .collect()
        .map { file ->
            [ [id:"domtblout"], file ]
        }

    ch_tsv = ch_families.tsv
        .map { meta, files ->
            files
        }
        .collect()
        .map { file ->
            [ [id:"tsv"], file ]
        }

    ch_discarded = ch_families.discarded
        .map { meta, files ->
            files
        }
        .collect()
        .map { file ->
            [ [id:"discarded"], file ]
        }

    ch_successful = ch_families.successful
        .map { meta, files ->
            files
        }
        .collect()
        .map { file ->
            [ [id:"successful"], file ]
        }

    ch_converged = ch_families.converged
        .map { meta, files ->
            files
        }
        .collect()
        .map { file ->
            [ [id:"converged"], file ]
        }

    ch_metadata = ch_families.metadata
        .map { meta, files ->
            files
        }
        .collect()
        .map { file ->
            [ [id:"metadata"], file ]
        }

    ch_logs = ch_families.logs
        .map { meta, files ->
            files
        }
        .collect()
        .map { file ->
            [ [id:"logs"], file ]
        }

    generated_families = REMOVE_REDUNDANCY( seed_msa_sto_dir, ch_seed_msa_sto, \
        ch_msa_sto, ch_hmm, \
        ch_rf, ch_domtblout, ch_tsv, \
        ch_discarded, ch_successful, ch_converged, \
        ch_metadata, ch_logs, ch_tm_ids, prob_ids, ch_rep_fa )
    ch_versions = ch_versions.mix( REMOVE_REDUNDANCY.out.versions )

    emit:
    versions     = ch_versions
    seed_msa_sto = generated_families.seed_msa_sto
    msa_sto      = generated_families.msa_sto
    metadata     = generated_families.metadata
    converged    = generated_families.converged
    tsv          = generated_families.tsv
    hmm          = generated_families.hmm
    rf           = generated_families.rf
}
