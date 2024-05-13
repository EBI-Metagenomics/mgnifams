#!/usr/bin/env nextflow

include { EXPORT_MGNIFAMS_CSV } from "${params.moduleDir}/export.nf"

workflow {
    EXPORT_MGNIFAMS_CSV(params.mgnifams_output_dir)
}
