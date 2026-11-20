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
include { EXPORT_FAMILIES_TSV            } from '../modules/local/export_families_tsv/main'
include { SPLIT_FAMILY_MEMBERS           } from '../modules/local/split_family_members/main'
include { BUILD_PARQUET_DOMAIN_QUERIES   } from '../modules/local/build_parquet_domain_queries/main'
include { PARSE_DOMAINS                  } from '../modules/local/parse_domains/main'
include { BIOMES_FROM_PARQUET            } from '../subworkflows/local/biomes_from_parquet'
include { afdbComputed                   } from '../subworkflows/local/utils_nfcore_mgnifams_pipeline'
include { batchQueryResults              } from '../subworkflows/local/utils_nfcore_mgnifams_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MGNIFAMS {

    take:
    ch_samplesheet // channel: [ meta, protein_input FASTA or [], sequences parquet or [], clusters parquet, pfam parquet, mgnifam_regions dir or [] ]
    parquet_chunks
    min_sequence_length
    outdir
    minimum_members
    clusters_chunk_size
    redundant_score_threshold
    redundant_residue_threshold
    starting_id
    pdb_chunk_size
    esmfold_db
    esmfold_params_path
    esmfold_3B_v1
    esm2_t36_3B_UR50D
    esm2_t36_3B_UR50D_contact_regression
    num_recycles_esmfold
    pdb_chunk_size_long
    skip_deeptmhmm
    deeptmhmm_path
    pfam_path
    funfams_path
    s4pred_chunk_size
    hh_mode
    hhdb_path
    foldseek_pdb_db
    foldseek_alphafold_db
    query_hmm_length_threshold
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    mgnifams_release        // string: this release, MAJOR.MINOR (families.tsv.gz)
    mgnify_proteins_release // string: MGnify Proteins release of the input, YYYY_MM
    family_query_chunks     // integer: about this many domain / biome lookup tasks, whole families each
    query_result_chunks     // integer: families per PARSE_DOMAINS / PARSE_BIOMES task

    main:

    ch_versions      = channel.empty()
    ch_multiqc_files = channel.empty()

    SETUP_CLUSTERS( ch_samplesheet, parquet_chunks, min_sequence_length, outdir, \
        minimum_members, clusters_chunk_size )
    ch_versions = ch_versions.mix( SETUP_CLUSTERS.out.versions )

    GENERATE_NONREDUNDANT_FAMILIES( SETUP_CLUSTERS.out.cluster_chunks, \
        SETUP_CLUSTERS.out.mgnifams_input_fa, redundant_score_threshold, redundant_residue_threshold, \
        starting_id )
    ch_versions = ch_versions.mix( GENERATE_NONREDUNDANT_FAMILIES.out.versions )

    // FTP families.tsv.gz: one `new` row per family (A17)
    EXPORT_FAMILIES_TSV(
        GENERATE_NONREDUNDANT_FAMILIES.out.metadata.map { meta, csv -> [ meta, [], [], csv ] },
        mgnifams_release,
        mgnify_proteins_release
    )
    ch_versions = ch_versions.mix( EXPORT_FAMILIES_TSV.out.versions )

    //
    // Members' domain architectures (Pfam parquet) and biomes (clusters parquet), as update_mgnifams
    // (SPLIT_FAMILY_MEMBERS fails if a family is not one contiguous run of rows)
    //
    SPLIT_FAMILY_MEMBERS( GENERATE_NONREDUNDANT_FAMILIES.out.tsv, family_query_chunks )
    ch_versions = ch_versions.mix( SPLIT_FAMILY_MEMBERS.out.versions )
    ch_family_chunks = SPLIT_FAMILY_MEMBERS.out.chunks
        .flatMap { _meta, chunks -> [chunks].flatten().collect { chunk -> [ [id: chunk.baseName], [ chunk ] ] } }

    BUILD_PARQUET_DOMAIN_QUERIES(
        ch_family_chunks.combine( ch_samplesheet.map { _meta, _fasta, _sequences, _clusters, pfam, _regions -> pfam } ),
        file(pfam_path, checkIfExists: true)
    )
    ch_versions = ch_versions.mix( BUILD_PARQUET_DOMAIN_QUERIES.out.versions )

    PARSE_DOMAINS(
        BUILD_PARQUET_DOMAIN_QUERIES.out.res.flatMap { meta, files -> batchQueryResults( meta, files, query_result_chunks ) },
        BUILD_PARQUET_DOMAIN_QUERIES.out.pfam_mapping.first(),
        GENERATE_NONREDUNDANT_FAMILIES.out.tsv.first()
    )
    ch_versions = ch_versions.mix( PARSE_DOMAINS.out.versions )

    BIOMES_FROM_PARQUET(
        ch_family_chunks,
        ch_samplesheet.map { _meta, _fasta, _sequences, clusters, _pfam, _regions -> clusters }.first(),
        query_result_chunks
    )
    ch_versions = ch_versions.mix( BIOMES_FROM_PARQUET.out.versions )

    PREDICT_STRUCTURES( GENERATE_NONREDUNDANT_FAMILIES.out.family_ids_fasta, \
        pdb_chunk_size, esmfold_db, esmfold_params_path, esmfold_3B_v1, \
        esm2_t36_3B_UR50D, esm2_t36_3B_UR50D_contact_regression, num_recycles_esmfold, \
        pdb_chunk_size_long, outdir )
    ch_versions = ch_versions.mix( PREDICT_STRUCTURES.out.versions )

    ANNOTATE_FAMILIES( GENERATE_NONREDUNDANT_FAMILIES.out.family_ids_fasta, \
        skip_deeptmhmm, deeptmhmm_path, pfam_path, funfams_path, s4pred_chunk_size, \
        GENERATE_NONREDUNDANT_FAMILIES.out.seed_msa, GENERATE_NONREDUNDANT_FAMILIES.out.full_msa, \
        hh_mode, hhdb_path, PREDICT_STRUCTURES.out.pdb, foldseek_pdb_db, foldseek_alphafold_db, outdir, starting_id )
    ch_versions = ch_versions.mix( ANNOTATE_FAMILIES.out.versions )

    //
    // What this run computed, for the DB build (same keys as update_mgnifams' update_info.csv)
    //
    def tm_computed = !skip_deeptmhmm && !workflow.profile.contains("conda")
    ANNOTATE_FAMILIES.out.z
        .map { z ->
            [
                'key,value',
                "tm_computed,${tm_computed}",
                "biome_computed,true",
                "afdb_computed,${afdbComputed()}",
                "pfam_lib,${file(pfam_path).name}",
                "funfams_lib,${file(funfams_path).name}",
                "mgnifams_release,${mgnifams_release}",
                "mgnify_proteins_release,${mgnify_proteins_release}",
                "pipeline_version,${workflow.manifest.version}",
                "hmmsearch_z,${z}",
                "delta_type,new_families"
            ].join('\n')
        }
        .collectFile( name: 'update_info.csv', storeDir: "${outdir}/generate_families", newLine: true )

    EXPORT_DATA( GENERATE_NONREDUNDANT_FAMILIES.out.metadata, PREDICT_STRUCTURES.out.scores, \
        ANNOTATE_FAMILIES.out.composition, ANNOTATE_FAMILIES.out.tm_composition, \
        ANNOTATE_FAMILIES.out.pfam_domains, ANNOTATE_FAMILIES.out.funfam_domains, query_hmm_length_threshold, \
        ANNOTATE_FAMILIES.out.pfam_model_hits, ANNOTATE_FAMILIES.out.foldseek_hits, outdir )
    ch_versions = ch_versions.mix( EXPORT_DATA.out.versions )

    //
    // Collate and save software versions
    //
    // nf-core modules emit [ process, tool, version ] tuples to the `versions` topic;
    // local modules still emit versions.yml files through ch_versions
    def topic_versions_string = channel.topic("versions")
        .distinct()
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            "${process}:\n${tool_versions.unique().sort().join('\n')}"
        }

    softwareVersionsToYAML(ch_versions)
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name: 'nf_core_pipeline_software_mqc_versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    summary_params      = paramsSummaryMap(parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = multiqc_methods_description ?
        file(multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    ch_multiqc_files = ch_multiqc_files.mix(SETUP_CLUSTERS.out.seqkit_stats_mqc.collect { t -> t[1] }.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(SETUP_CLUSTERS.out.cluster_distr_mqc.collect { t -> t[1] }.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(GENERATE_NONREDUNDANT_FAMILIES.out.discarded_mqc.collect { t -> t[1] }.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(GENERATE_NONREDUNDANT_FAMILIES.out.metadata_mqc.collect { t -> t[1] }.ifEmpty([]))

    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'mgnifams'],
                files,
                // The custom config goes last, so its settings take precedence over the default one
                [ file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true) ] +
                    (multiqc_config ? [ file(multiqc_config, checkIfExists: true) ] : []),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )

    emit:
    multiqc_report = MULTIQC.out.report.map { _meta, report -> report }.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
