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
include { EXPORT_DB                      } from "../subworkflows/export_db.nf"

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
    
    setup_res          = SETUP_CLUSTERS(ch_input, params.fasta_input_mode, params.compress_mode)    
    generated_families = GENERATE_NONREDUNDANT_FAMILIES(setup_res.clusters_tsv, [], setup_res.mgnifams_input_fa)
    annotated_families = ANNOTATE_FAMILIES(generated_families.seed_msa_sto, generated_families.msa_sto)
    EXPORT_DB(generated_families.metadata, generated_families.converged, generated_families.tsv, \
        annotated_families.pfam_hits, annotated_families.foldseek_hits, annotated_families.scores_ch, \
        annotated_families.cif_ch, annotated_families.fa_seed_msa_ch, annotated_families.fa_msa_ch, \
        generated_families.hmm, generated_families.rf, params.fasta_input_mode)

    CUSTOM_DUMPSOFTWAREVERSIONS(channel.topic('versions').unique().collectFile(name: 'collated_versions.yml'))

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/