#!/usr/bin/env nextflow

include { annotate_structures } from "$baseDir/subworkflows/annotate_structures/main.nf"

workflow {
    Channel
        .fromPath(params.wp1_unannotated_ids_path) 
        .set { wp1_unannotated_ids }

    Channel
        .fromPath(params.family_reps_path) 
        .set { family_reps }

    annotate_structures(wp1_unannotated_ids, family_reps)
    
}