#!/usr/bin/env nextflow

include { GENERATE_FAMILIES_PARALLEL     } from "${projectDir}/../../subworkflows/generate_families_parallel/main.nf"
include { MOVE_TO_DIR as MOVE_SEED_MSA   } from "${projectDir}/../../modules/family/main.nf"
include { MOVE_TO_DIR as MOVE_MSA        } from "${projectDir}/../../modules/family/main.nf"
include { MOVE_TO_DIR as MOVE_HMM        } from "${projectDir}/../../modules/family/main.nf"
include { MOVE_TO_DIR as MOVE_RF         } from "${projectDir}/../../modules/family/main.nf"
include { MOVE_TO_DIR as MOVE_DOMTBL     } from "${projectDir}/../../modules/family/main.nf"
include { MOVE_TO_DIR as MOVE_TSV        } from "${projectDir}/../../modules/family/main.nf"
include { MOVE_TO_DIR as MOVE_DISCARDED  } from "${projectDir}/../../modules/family/main.nf"
include { MOVE_TO_DIR as MOVE_SUCCESSFUL } from "${projectDir}/../../modules/family/main.nf"
include { MOVE_TO_DIR as MOVE_CONVERGED  } from "${projectDir}/../../modules/family/main.nf"
include { MOVE_TO_DIR as MOVE_METADATA   } from "${projectDir}/../../modules/family/main.nf"
include { MOVE_TO_DIR as MOVE_LOGS       } from "${projectDir}/../../modules/family/main.nf"
include { MOVE_TO_DIR as MOVE_FAMILY     } from "${projectDir}/../../modules/family/main.nf"
include { REMOVE_REDUNDANCY              } from "${projectDir}/../../subworkflows/remove_redundancy/main.nf"

workflow {
    families_ch = GENERATE_FAMILIES_PARALLEL( params.clusters_tsv, params.checked_clusters_txt, params.mgnifams_input_fasta_path)

    seed_msa_sto_ch = families_ch.seed_msa_sto.collect()
    msa_sto_ch      = families_ch.msa_sto.collect()
    hmm_ch          = families_ch.hmm.collect()
    rf_ch           = families_ch.rf.collect()
    domtblout_ch    = families_ch.domtblout.collect()
    tsv_ch          = families_ch.tsv.collect()
    discarded_ch    = families_ch.discarded.collect()
    successful_ch   = families_ch.successful.collect()
    converged_ch    = families_ch.converged.collect()
    metadata_ch     = families_ch.metadata.collect()
    logs_ch         = families_ch.logs.collect()

    seed_msa_sto_dir = MOVE_SEED_MSA(seed_msa_sto_ch, "seed_msa_sto")
    msa_sto_dir      = MOVE_MSA(msa_sto_ch, "msa_sto")
    hmm_dir          = MOVE_HMM(hmm_ch, "hmm")
    rf_dir           = MOVE_RF(rf_ch, "rf")
    domtblout_dir    = MOVE_DOMTBL(domtblout_ch, "domtblout")
    tsv_dir          = MOVE_TSV(tsv_ch, "refined_families")
    discarded_dir    = MOVE_DISCARDED(discarded_ch, "discarded_clusters")
    successful_dir   = MOVE_SUCCESSFUL(successful_ch, "successful_clusters")
    converged_dir    = MOVE_CONVERGED(converged_ch, "converged_families")
    metadata_dir     = MOVE_METADATA(metadata_ch, "family_metadata")
    logs_dir         = MOVE_LOGS(logs_ch, "logs")

    families_all_dir = MOVE_FAMILY(seed_msa_sto_dir.mix(msa_sto_dir, \
        hmm_dir, rf_dir, domtblout_dir, tsv_dir, \
        discarded_dir, successful_dir, converged_dir, \
        metadata_dir, logs_dir).collect(), "families_prepooled")

    seed_msa_sto_dir
        .map { filepath ->
            return [ [id:"remove_redundancy"], file(filepath) ]
        }
        .set { seed_msa_sto_ch }
    REMOVE_REDUNDANCY(seed_msa_sto_ch, families_all_dir)
}
