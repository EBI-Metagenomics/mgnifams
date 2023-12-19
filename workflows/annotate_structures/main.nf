#!/usr/bin/env nextflow

include { ANNOTATE_STRUCTURES } from "$launchDir/subworkflows/annotate_structures/main.nf"

workflow {
    pdb = [ [ id:'pdb_data' ], params.pdb_path ]
    ANNOTATE_STRUCTURES(pdb)
}
