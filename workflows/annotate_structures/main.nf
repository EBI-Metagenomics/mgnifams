#!/usr/bin/env nextflow

include { ANNOTATE_STRUCTURES } from "${projectDir}/../../subworkflows/annotate_structures/main.nf"

workflow {
    pdb_ch = [ [ id:'pdb_data' ], params.pdb_path ]
    ANNOTATE_STRUCTURES(pdb_ch)
}
