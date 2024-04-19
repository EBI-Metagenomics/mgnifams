#!/usr/bin/env nextflow

include { EXPORT_MGNIFAMS_CSV } from "$launchDir/modules/export.nf"

workflow EXPORT_TABLES {
    take:
    mgnifams_output_dir

    main:
    EXPORT_MGNIFAMS_CSV(mgnifams_output_dir)
}
