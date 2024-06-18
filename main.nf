#!/usr/bin/env nextflow

include { validateParameters; paramsHelp; paramsSummaryLog } from 'plugin/nf-schema'

if (params.help) {
    log.info paramsHelp("nextflow run main.nf -c conf/end-to-end.config -profile slurm")
    exit 0
}
validateParameters()
log.info paramsSummaryLog(workflow)

include { PREPROCESS_INPUT                      } from "${projectDir}/subworkflows/preprocess_input/main.nf"
include { INITIATE_PROTEINS                     } from "${projectDir}/subworkflows/initiate_proteins/main.nf"
include { EXECUTE_CLUSTERING                    } from "${projectDir}/subworkflows/execute_clustering/main.nf"
include { GENERATE_FAMILIES_PARALLEL            } from "${projectDir}/subworkflows/generate_families_parallel/main.nf"
include { MOVE_TO_DIR                           } from "${projectDir}/modules/family/main.nf"
include { REMOVE_REDUNDANCY                     } from "${projectDir}/subworkflows/remove_redundancy/main.nf"
include { REFORMAT_MSA as REFORMAT_SEED_MSA     } from "${projectDir}/subworkflows/reformat_msa/main.nf"
include { REFORMAT_MSA as REFORMAT_HMMALIGN_MSA } from "${projectDir}/subworkflows/reformat_msa/main.nf"
include { ANNOTATE_MODELS                       } from "${projectDir}/subworkflows/annotate_models/main.nf"
include { PREDICT_STRUCTURES                    } from "${projectDir}/subworkflows/predict_structures/main.nf"
include { ANNOTATE_STRUCTURES                   } from "${projectDir}/subworkflows/annotate_structures/main.nf"

workflow {
    preprocessed_sequence_explorer_protein_ch = PREPROCESS_INPUT(params.sequence_explorer_protein_path, params.compress_mode).preprocessed_sequence_explorer_protein_ch
    fasta_ch                                  = INITIATE_PROTEINS( preprocessed_sequence_explorer_protein_ch ).fasta_ch
    clusters                                  = EXECUTE_CLUSTERING( fasta_ch )
    clusters_tsv                              = clusters.clusters_tsv.map { meta, filepath -> filepath }
    starting_num_sequences                    = clusters.num_sequences
    // TODO update below here
    generated_families                        = GENERATE_FAMILIES_ALL(clusters_tsv, fasta_ch, starting_num_sequences)

    generated_families.seed_msa_sto
        .map { files ->
            String filePath = files[0]
            int lastIndex = filePath.lastIndexOf('/')
            String seed_msa_dir = filePath.substring(0, lastIndex + 1)
            return [ [id:"seed_msa"], file(seed_msa_dir) ]
        }
        .set { seed_msa_ch }
    
    generated_families.msa_sto
        .map { files ->
            String filePath = files[0]
            int lastIndex = filePath.lastIndexOf('/')
            String msa_dir = filePath.substring(0, lastIndex + 1)
            return [ [id:"msa"], file(msa_dir) ]
        }
        .set { hmmalign_msa_ch }

    fa_seed_msa_ch = REFORMAT_SEED_MSA(seed_msa_ch).fa_ch
    REFORMAT_HMMALIGN_MSA( hmmalign_msa_ch )
    ANNOTATE_MODELS( fa_seed_msa_ch )
    pdb_ch = PREDICT_STRUCTURES(hmmalign_msa_ch).pdb_ch
    ANNOTATE_STRUCTURES(pdb_ch)
}
