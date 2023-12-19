#!/usr/bin/env nextflow

include { ANNOTATE_STRUCTURES } from "$launchDir/subworkflows/annotate_structures/main.nf"

workflow {
    ANNOTATE_STRUCTURES(params.pdb_path)
}
