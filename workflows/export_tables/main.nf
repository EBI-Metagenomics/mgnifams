#!/usr/bin/env nextflow

include { EXPORT_TABLES } from "$launchDir/subworkflows/export_tables/main.nf"

workflow {
    EXPORT_TABLES(params.mgnifams_output_dir)
}
