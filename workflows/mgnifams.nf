/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_mgnifams_pipeline'

include { SETUP_CLUSTERS                 } from '../subworkflows/local/setup_clusters'
include { GENERATE_NONREDUNDANT_FAMILIES } from '../subworkflows/local/generate_nonredundant_families'
include { PREDICT_STRUCTURES             } from '../subworkflows/local/predict_structures'
include { ANNOTATE_FAMILIES              } from '../subworkflows/local/annotate_families'
include { EXPORT_DATA                    } from '../subworkflows/local/export_data'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MGNIFAMS {
    
    take:
    ch_samplesheet // channel: samplesheet read in from --input
    fasta_input_mode
    compress_mode
    input_csv_chunk_size
    min_sequence_length
    outdir
    minimum_members
    clusters_chunk_size
    mgnifams_discard_min_rep_length
    mgnifams_discard_max_rep_length
    mgnifams_discard_min_starting_membership
    mgnifams_max_seq_identity
    mgnifams_max_seed_seqs
    mgnifams_max_gap_occupancy
    mgnifams_recruit_evalue_cutoff
    mgnifams_recruit_hit_length_percentage
    redundant_length_threshold
    redundant_score_threshold
    similarity_score_threshold
    starting_id
    fold_mode
    pdb_chunk_size
    esmfold_db
    esmfold_params_path
    esmfold_3B_v1
    esm2_t36_3B_UR50D
    esm2_t36_3B_UR50D_contact_regression
    num_recycles_esmfold
    pdb_chunk_size_long
    alphafold3_db
    alphafold3_path
    alphafold3_small_bfd_path
    alphafold3_mgnify_path
    alphafold3_pdb_mmcif_path
    alphafold3_uniref90_path
    alphafold3_pdb_seqres_path
    alphafold3_uniprot_path
    alphafold3_small_bfd_link
    alphafold3_mgnify_link
    alphafold3_pdb_mmcif_link
    alphafold3_uniref90_link
    alphafold3_pdb_seqres_link
    uniprot_link
    funfams_path
    hh_mode
    hhdb_path
    foldseek_db_path
    multiqc_config
    multiqc_logo
    multiqc_methods_description

    main:
    
    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()

    SETUP_CLUSTERS( ch_samplesheet, fasta_input_mode, compress_mode, \
        input_csv_chunk_size, min_sequence_length, outdir, \
        minimum_members, clusters_chunk_size )
    ch_versions = ch_versions.mix( SETUP_CLUSTERS.out.versions )

    GENERATE_NONREDUNDANT_FAMILIES( SETUP_CLUSTERS.out.cluster_chunks, \
        SETUP_CLUSTERS.out.mgnifams_input_fa,  mgnifams_discard_min_rep_length, \
        mgnifams_discard_max_rep_length, mgnifams_discard_min_starting_membership, \
        mgnifams_max_seq_identity, mgnifams_max_seed_seqs, mgnifams_max_gap_occupancy, \
        mgnifams_recruit_evalue_cutoff, mgnifams_recruit_hit_length_percentage, \
        outdir, redundant_length_threshold, redundant_score_threshold, \
        similarity_score_threshold, starting_id )
    ch_versions = ch_versions.mix( GENERATE_NONREDUNDANT_FAMILIES.out.versions )

    PREDICT_STRUCTURES( GENERATE_NONREDUNDANT_FAMILIES.out.family_reps, fold_mode, \
        pdb_chunk_size, esmfold_db, esmfold_params_path, esmfold_3B_v1, \
        esm2_t36_3B_UR50D, esm2_t36_3B_UR50D_contact_regression, num_recycles_esmfold, \
        pdb_chunk_size_long, outdir, alphafold3_db, alphafold3_path, alphafold3_small_bfd_path, \
        alphafold3_mgnify_path, alphafold3_pdb_mmcif_path, alphafold3_uniref90_path,\
        alphafold3_pdb_seqres_path, alphafold3_uniprot_path, alphafold3_small_bfd_link,\
        alphafold3_mgnify_link, alphafold3_pdb_mmcif_link, alphafold3_uniref90_link,\
        alphafold3_pdb_seqres_link, uniprot_link )
    ch_versions = ch_versions.mix( PREDICT_STRUCTURES.out.versions )
    
    ANNOTATE_FAMILIES( GENERATE_NONREDUNDANT_FAMILIES.out.family_reps, funfams_path, \
        GENERATE_NONREDUNDANT_FAMILIES.out.seed_msa, GENERATE_NONREDUNDANT_FAMILIES.out.full_msa, \
        hh_mode, hhdb_path, PREDICT_STRUCTURES.out.pdb, foldseek_db_path, outdir )
    ch_versions = ch_versions.mix( ANNOTATE_FAMILIES.out.versions )

    // TODO update vs remove
    // EXPORT_DATA( GENERATE_NONREDUNDANT_FAMILIES.out.metadata, GENERATE_NONREDUNDANT_FAMILIES.out.converged, \
    //    GENERATE_NONREDUNDANT_FAMILIES.out.tsv, ANNOTATE_FAMILIES.out.pfam_hits, \
    //    ANNOTATE_FAMILIES.out.foldseek_hits, ANNOTATE_FAMILIES.out.scores )
    // ch_versions = ch_versions.mix( EXPORT_DATA.out.versions )

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name: 'nf_core_pipeline_software_mqc_versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = multiqc_config ?
        Channel.fromPath(multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = multiqc_logo ?
        Channel.fromPath(multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = multiqc_methods_description ?
        file(multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    ch_multiqc_files = ch_multiqc_files.mix(SETUP_CLUSTERS.out.seqkit_stats_mqc.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(SETUP_CLUSTERS.out.cluster_distr_mqc.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(GENERATE_NONREDUNDANT_FAMILIES.out.discarded_mqc.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(GENERATE_NONREDUNDANT_FAMILIES.out.metadata_mqc.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(GENERATE_NONREDUNDANT_FAMILIES.out.similarity_mqc.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/