#!/usr/bin/env nextflow

include { CREATE_CLUSTERS_PKL } from "${launchDir}/modules/family/main.nf"
include { CREATE_FASTA_PKL    } from "${launchDir}/modules/family/main.nf"
include { REFINE_FAMILIES     } from "${launchDir}/modules/family/main.nf"

workflow GENERATE_FAMILIES {
    take:
    families_tsv
    mgnifams_fasta
    mode
    
    main:
    if (mode == "load_state") {
        families_pkl = file(params.families_pkl_path)
        mgnifams_fa = file(params.mgnifams_input_fasta_path)
        mgnifams_pkl = file(params.mgnifams_input_pkl_path)
    } else {
        families_pkl = CREATE_CLUSTERS_PKL(families_tsv).pkl
        mgnifams_fa = mgnifams_fasta
        mgnifams_pkl = CREATE_FASTA_PKL(mgnifams_fasta).pkl
    }
    refined_families = REFINE_FAMILIES(families_pkl, mgnifams_fa, mgnifams_pkl, params.minimum_members)

    emit:
    tsv = refined_families.tsv
    seed_msa = refined_families.seed_msa
    msa = refined_families.msa
    hmm = refined_families.hmm
    domtblout = refined_families.domtblout
}
