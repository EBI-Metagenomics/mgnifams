/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog } from 'plugin/nf-schema'

log.info paramsSummaryLog(workflow)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// setup_clusters
include { PREPROCESS_INPUT   } from "../subworkflows/preprocess_input/main.nf"
include { INITIATE_PROTEINS  } from "../subworkflows/initiate_proteins/main.nf"
include { EXECUTE_CLUSTERING } from "../subworkflows/execute_clustering/main.nf"

// generate_nonredundant_families
include { GENERATE_FAMILIES_PARALLEL } from "../subworkflows/generate_families_parallel/main.nf"
include { FLAG_TRANSMEMBRANE         } from "../subworkflows/flag_transmembrane/main.nf"
include { MOVE_TO_DIR                } from "../modules/family/main.nf"
include { REMOVE_REDUNDANCY          } from "../subworkflows/remove_redundancy/main.nf"

// annotate_families
include { REFORMAT_MSA as REFORMAT_SEED_MSA     } from "../subworkflows/reformat_msa/main.nf"
include { REFORMAT_MSA as REFORMAT_HMMALIGN_MSA } from "../subworkflows/reformat_msa/main.nf"
include { ANNOTATE_MODELS                       } from "../subworkflows/annotate_models/main.nf"
include { PREDICT_STRUCTURES                    } from "../subworkflows/predict_structures/main.nf"
include { ANNOTATE_STRUCTURES                   } from "../subworkflows/annotate_structures/main.nf"

// export_db
include { EXPORT_DB } from "../subworkflows/export_db/main.nf"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MGNIFAMS {
    // setup_clusters
    processed_input_protein_ch = PREPROCESS_INPUT(params.sequence_explorer_protein_path, params.compress_mode).processed_input_protein_ch
    mgnifams_input_fasta_ch    = INITIATE_PROTEINS( processed_input_protein_ch ).fasta_ch
    clusters                   = EXECUTE_CLUSTERING( mgnifams_input_fasta_ch )
    clusters_tsv               = clusters.clusters_tsv.map { meta, filepath -> filepath }

    // generate_nonredundant_families
    families_ch = GENERATE_FAMILIES_PARALLEL( clusters_tsv, [], mgnifams_input_fasta_ch)

    msa_sto_ch = families_ch.msa_sto.collect()
    msa_sto_ch
        .map { files ->
            return [ [id:"flag_tm"], files ]
        }
        .set { tm_msa_ch }
    
    tm_ch     = FLAG_TRANSMEMBRANE(tm_msa_ch)
    rep_fa_ch = tm_ch.fa_ch
    tm_ids_ch = tm_ch.tm_ids_ch
    
    seed_msa_sto_ch = families_ch.seed_msa_sto.collect()
    seed_msa_sto_dir = MOVE_TO_DIR(seed_msa_sto_ch, "seed_msa_sto")
    seed_msa_sto_dir
        .map { filepath ->
            return [ [id:"remove_redundancy"], file(filepath) ]
        }
        .set { seed_msa_sto_dir }
        
    generated_families = REMOVE_REDUNDANCY(seed_msa_sto_dir, seed_msa_sto_ch, \
        msa_sto_ch, families_ch.hmm.collect(), \
        families_ch.rf.collect(), families_ch.domtblout.collect(), families_ch.tsv.collect(), \
        families_ch.discarded.collect(), families_ch.successful.collect(), families_ch.converged.collect(), \
        families_ch.metadata.collect(), families_ch.logs.collect(), tm_ids_ch, rep_fa_ch)

    // annotate_families
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
    fa_msa_ch      = REFORMAT_HMMALIGN_MSA(hmmalign_msa_ch).fa_ch
    pfam_hits      = ANNOTATE_MODELS(fa_seed_msa_ch)
    structure_ch   = PREDICT_STRUCTURES(hmmalign_msa_ch)
    pdb_ch         = structure_ch.pdb_ch
    scores_ch      = structure_ch.scores_ch
    cif_ch         = structure_ch.cif_ch
    foldseek_hits  = ANNOTATE_STRUCTURES(pdb_ch)

    // export_db
    EXPORT_DB(generated_families.metadata, generated_families.converged, \
        generated_families.tsv, pfam_hits, foldseek_hits, scores_ch, \
        cif_ch, fa_seed_msa_ch, fa_msa_ch, \
        generated_families.hmm, generated_families.rf)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/