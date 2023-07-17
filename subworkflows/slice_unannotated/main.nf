#!/usr/bin/env nextflow

include { SLICE; SCANSLICE } from "$baseDir/modules/slicing/slicing.nf"

workflow slice_unannotated {
    take:
    msa_ch
    tblout_ch
    
    main:
    slice_ch = SLICE(tblout_ch)
    SCANSLICE(slice_ch, msa_ch)

    // emit:

}