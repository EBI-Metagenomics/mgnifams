#!/usr/bin/env nextflow

include { EXPORT_KNOWN_ANNOTATIONS_CSV; CONCAT_KNOWN_ANNOTATIONS } from "$baseDir/modules/exporting.nf"

workflow annotate_known {
    take:
    known_fasta
    
    main:
    known_annotations_ch = EXPORT_KNOWN_ANNOTATIONS_CSV(known_fasta)
    CONCAT_KNOWN_ANNOTATIONS(known_annotations_ch.collect())
}