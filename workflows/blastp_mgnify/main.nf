#!/usr/bin/env nextflow

include { BLASTP_MGNIFY } from "$launchDir/subworkflows/blastp_mgnify/main.nf"

workflow {
    BLASTP_MGNIFY(params.anti_defence_proteins_path, params.mgnifams_input_path, params.clusters_tsv)
}
