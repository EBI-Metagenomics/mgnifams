#!/usr/bin/env nextflow

include { PREPROCESS_INPUT    } from "$launchDir/subworkflows/preprocess_input/main.nf"
include { INITIATE_PROTEINS   } from "$launchDir/subworkflows/initiate_proteins/main.nf"
include { EXECUTE_CLUSTERING  } from "$launchDir/subworkflows/execute_clustering/main.nf"
include { GENERATE_FAMILIES   } from "$launchDir/subworkflows/generate_families/main.nf"
include { ANNOTATE_MODELS     } from "$launchDir/subworkflows/annotate_models/main.nf"
include { PREDICT_STRUCTURES  } from "$launchDir/subworkflows/predict_structures/main.nf"
include { ANNOTATE_STRUCTURES } from "$launchDir/subworkflows/annotate_structures/main.nf"

workflow {
    Channel
        .fromPath(params.mgy90_path)
        .set { mgy90_file_bz2 }

    preprocessed_mgy90_file = PREPROCESS_INPUT( mgy90_file_bz2 ).preprocessed_mgy90_file
    out_fasta               = INITIATE_PROTEINS( preprocessed_mgy90_file ).out_fasta
    families_tsv            = EXECUTE_CLUSTERING( out_fasta ).families_tsv.map { meta, families_tsv -> families_tsv }
    generated_families      = GENERATE_FAMILIES( families_tsv, params.empty_file, out_fasta )
    generated_families.seed_msa
        .map { files ->
            String filePath = files[0]
            int lastIndex = filePath.lastIndexOf('/')
            String seed_msa_dir = filePath.substring(0, lastIndex + 1)
            return [ [id:"annotated_models"], file(seed_msa_dir) ]
        }
        .set { seed_msa_dir }
    unannotated = ANNOTATE_MODELS( seed_msa_dir, params.hhdb_folder_path, "hhblits" )
    generated_families.msa
        .map { files ->
            String filePath = files[0]
            int lastIndex = filePath.lastIndexOf('/')
            String msa_dir = filePath.substring(0, lastIndex + 1)
        }
        .set { msa_dir }
    pdb_ch = PREDICT_STRUCTURES( msa_dir, unannotated, "all" ).pdb_ch
    ANNOTATE_STRUCTURES( pdb_ch )
}
