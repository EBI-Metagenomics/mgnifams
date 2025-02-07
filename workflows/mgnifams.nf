/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog } from 'plugin/nf-schema'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SETUP_CLUSTERS                 } from "../subworkflows/setup_clusters.nf"
include { GENERATE_NONREDUNDANT_FAMILIES } from "../subworkflows/generate_nonredundant_families.nf"
include { ANNOTATE_FAMILIES              } from "../subworkflows/annotate_families.nf"
include { EXPORT_DATA                    } from "../subworkflows/export_data.nf"

include { CUSTOM_DUMPSOFTWAREVERSIONS    } from '../modules/nf-core/custom/dumpsoftwareversions/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

nextflow.preview.topic = true

workflow MGNIFAMS {
    
    take:
    ch_input

    main:
    log.info paramsSummaryLog(workflow)

    SETUP_CLUSTERS(ch_input, params.fasta_input_mode, params.compress_mode)
    generated_families = GENERATE_NONREDUNDANT_FAMILIES(SETUP_CLUSTERS.out.clusters_tsv, [], SETUP_CLUSTERS.out.mgnifams_input_fa)
    annotated_families = ANNOTATE_FAMILIES(generated_families.seed_msa_sto, generated_families.msa_sto)
    EXPORT_DATA(generated_families.metadata, generated_families.converged, generated_families.tsv, \
        annotated_families.pfam_hits, annotated_families.foldseek_hits, annotated_families.scores_ch)

    CUSTOM_DUMPSOFTWAREVERSIONS(channel.topic('versions').unique().collectFile(name: 'collated_versions.yml'))

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/