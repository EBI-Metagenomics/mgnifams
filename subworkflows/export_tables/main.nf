#!/usr/bin/env nextflow

include { EXPORT_PROTEINS_CSV } from "$launchDir/modules/export.nf"
include { EXPORT_MGNIFAMS_CSV } from "$launchDir/modules/export.nf"

workflow EXPORT_TABLES {
    take:
    mgnifams_input_path
    mgnifams_output_dir

    main:
    EXPORT_PROTEINS_CSV(mgnifams_input_path)
    print(mgnifams_output_dir)
    EXPORT_MGNIFAMS_CSV(mgnifams_output_dir)
}
