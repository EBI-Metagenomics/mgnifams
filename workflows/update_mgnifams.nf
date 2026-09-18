/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    UPDATE_MGNIFAMS: refresh full MSAs and representatives of existing families
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Seed MSAs and HMMs are not touched (update_families --skip_refine). Everything derived
    from the family representatives is recomputed and loaded into a delta sqlite DB with
    only the successful families, ready to be merged (assets/merge_update_delta.sql).
----------------------------------------------------------------------------------------
*/

include { MULTIQC                      } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap             } from 'plugin/nf-schema'
include { paramsSummaryMultiqc         } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML       } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText       } from '../subworkflows/local/utils_nfcore_mgnifams_pipeline'

include { UPDATE_FAMILIES              } from '../subworkflows/local/update_families'
include { PREDICT_STRUCTURES           } from '../subworkflows/local/predict_structures'
include { ANNOTATE_REPS                } from '../subworkflows/local/annotate_reps'
include { ANNOTATE_STRUCTURES          } from '../subworkflows/local/annotate_structures'
include { EXPORT_DATA                  } from '../subworkflows/local/export_data'
include { BUILD_PARQUET_DOMAIN_QUERIES } from '../modules/local/build_parquet_domain_queries/main'
include { PARSE_DOMAINS                } from '../modules/local/parse_domains/main'
include { QUERY_MGNPROTEIN_DB          } from '../modules/local/query_mgnprotein_db/main'
include { PARSE_BIOMES                 } from '../modules/local/parse_biomes/main'
include { INIT_SQLITE                  } from '../modules/local/init_sqlite/main'
include { IMPORT_QUERIES               } from '../modules/local/import_queries/main'
include { UPDATE_SQLITE_BLOBS_STAGED   } from '../modules/local/update_sqlite_blobs_staged/main'
include { HHSUITE_REFORMAT as REFORMAT_FULL_MSA_A3M } from '../modules/nf-core/hhsuite/reformat/main'
include { TRUNCATE_A3M                 } from '../modules/local/truncate_a3m/main'
include { COLABFOLD_BATCH_MSA          } from '../modules/local/colabfold_batch_msa/main'

workflow UPDATE_MGNIFAMS {

    take:
    ch_samplesheet      // channel: [ meta, sequences_parquet, pfam_parquet, hmm_lib, db_config or [] ]
    parquet_chunks
    min_sequence_length
    hmm_chunk_size
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
    foldseek_db_path
    query_hmm_length_threshold
    query_result_chunks
    run_alphafold2
    colabfold_params_path
    af2_max_msa_seqs
    af2_num_recycles
    outdir
    multiqc_config
    multiqc_logo
    multiqc_methods_description

    main:
    ch_versions      = channel.empty()
    ch_multiqc_files = channel.empty()
    ch_meta          = ch_samplesheet.map { meta, _sequences, _pfam, _hmms, _db_config -> meta }
    // Same condition ANNOTATE_REPS uses to run DeepTMHMM
    def tm_computed  = !skip_deeptmhmm && !workflow.profile.contains("conda")

    UPDATE_FAMILIES( ch_samplesheet, parquet_chunks, min_sequence_length, hmm_chunk_size, outdir )
    ch_versions = ch_versions.mix( UPDATE_FAMILIES.out.versions )

    PREDICT_STRUCTURES( UPDATE_FAMILIES.out.family_ids_fasta, pdb_chunk_size, esmfold_db, esmfold_params_path, \
        esmfold_3B_v1, esm2_t36_3B_UR50D, esm2_t36_3B_UR50D_contact_regression, num_recycles_esmfold, \
        pdb_chunk_size_long, outdir )
    ch_versions = ch_versions.mix( PREDICT_STRUCTURES.out.versions )

    // Not ANNOTATE_FAMILIES: its HH-suite model annotation works on the seed MSA, which does not change
    //
    // Optional: AlphaFold2 (ColabFold) from each family's full MSA; published only, not used downstream
    //
    if (run_alphafold2) {
        REFORMAT_FULL_MSA_A3M( UPDATE_FAMILIES.out.full_msa, 'sto', 'a3m' )
        ch_versions = ch_versions.mix( REFORMAT_FULL_MSA_A3M.out.versions )

        TRUNCATE_A3M( REFORMAT_FULL_MSA_A3M.out.msa, af2_max_msa_seqs )
        ch_versions = ch_versions.mix( TRUNCATE_A3M.out.versions )

        ch_af2_batches = TRUNCATE_A3M.out.a3m
            .map { _meta, a3m -> a3m }
            .collect()
            .flatMap { a3ms ->
                a3ms.sort { f -> f.name }.collate( pdb_chunk_size ).withIndex().collect { batch, index -> [ [id: "af2_batch_${index}"], batch ] }
            }
        COLABFOLD_BATCH_MSA( ch_af2_batches, colabfold_params_path ? file(colabfold_params_path, checkIfExists: true) : [], af2_num_recycles )
        ch_versions = ch_versions.mix( COLABFOLD_BATCH_MSA.out.versions )
    }

    ANNOTATE_REPS( UPDATE_FAMILIES.out.family_ids_fasta, skip_deeptmhmm, deeptmhmm_path, pfam_path, funfams_path )
    ch_versions = ch_versions.mix( ANNOTATE_REPS.out.versions )

    ANNOTATE_STRUCTURES( PREDICT_STRUCTURES.out.pdb, foldseek_db_path, outdir )
    ch_versions = ch_versions.mix( ANNOTATE_STRUCTURES.out.versions )

    EXPORT_DATA( UPDATE_FAMILIES.out.metadata, PREDICT_STRUCTURES.out.scores, \
        ANNOTATE_REPS.out.composition, ANNOTATE_REPS.out.tm_composition, channel.value([ [id: 'seed_sizes'], [] ]), \
        ANNOTATE_REPS.out.pfam_domains, ANNOTATE_REPS.out.funfam_domains, query_hmm_length_threshold, \
        channel.empty(), ANNOTATE_STRUCTURES.out.foldseek_hits, outdir )
    ch_versions = ch_versions.mix( EXPORT_DATA.out.versions )

    //
    // Domain architectures of the new members, from the Pfam parquet
    //
    BUILD_PARQUET_DOMAIN_QUERIES(
        UPDATE_FAMILIES.out.refined_families.combine( ch_samplesheet.map { _meta, _sequences, pfam, _hmms, _db_config -> pfam } ),
        file(pfam_path, checkIfExists: true)
    )
    ch_versions = ch_versions.mix( BUILD_PARQUET_DOMAIN_QUERIES.out.versions )

    ch_domain_batches = BUILD_PARQUET_DOMAIN_QUERIES.out.res
        .flatMap { _meta, files ->
            [files].flatten().collate( query_result_chunks ).withIndex().collect { flist, index -> [ [id: "batch_${index}"], flist ] }
        }
    PARSE_DOMAINS( ch_domain_batches, BUILD_PARQUET_DOMAIN_QUERIES.out.pfam_mapping.first(), UPDATE_FAMILIES.out.refined_families.first() )
    ch_versions = ch_versions.mix( PARSE_DOMAINS.out.versions )

    //
    // Biome distributions of the new members, from the MGnify proteins DB (only with mgnprotein_db_config)
    //
    ch_db_config = ch_samplesheet
        .filter { _meta, _sequences, _pfam, _hmms, db_config -> db_config }
        .map { meta, _sequences, _pfam, _hmms, db_config -> [ meta, db_config ] }
    QUERY_MGNPROTEIN_DB( ch_db_config.combine( UPDATE_FAMILIES.out.refined_families.map { _meta, tsv -> tsv } ) )
    ch_versions = ch_versions.mix( QUERY_MGNPROTEIN_DB.out.versions )

    ch_biome_batches = QUERY_MGNPROTEIN_DB.out.res
        .flatMap { _meta, files ->
            [files].flatten().collate( query_result_chunks ).withIndex().collect { flist, index -> [ [id: "batch_${index}"], flist ] }
        }
    PARSE_BIOMES( ch_biome_batches, QUERY_MGNPROTEIN_DB.out.biome_mapping.first() )
    ch_versions = ch_versions.mix( PARSE_BIOMES.out.versions )

    //
    // Delta DB: tables of the successful families, then their blobs
    //
    INIT_SQLITE( ch_meta.map { meta -> [ meta, file("${projectDir}/assets/data/db_schema.sqlite", checkIfExists: true) ] } )
    ch_versions = ch_versions.mix( INIT_SQLITE.out.versions )

    ch_tables = EXPORT_DATA.out.mgnifam
        .mix( EXPORT_DATA.out.pfams, EXPORT_DATA.out.funfams, EXPORT_DATA.out.folds )
        .map { _meta, csv -> csv }
        .collect()
    IMPORT_QUERIES( ch_meta.combine( ch_tables.map { csvs -> [ csvs ] } ), INIT_SQLITE.out.db )
    ch_versions = ch_versions.mix( IMPORT_QUERIES.out.versions )

    ch_update_info = ch_samplesheet.map { _meta, _sequences, _pfam, _hmms, db_config ->
        [
            tm_computed     : tm_computed,
            biome_computed  : db_config ? true : false,
            child_tables    : 'pfams,funfams,folds',
            pfam_lib        : file(pfam_path).name,
            funfams_lib     : file(funfams_path).name,
            pipeline_version: workflow.manifest.version
        ]
    }
    UPDATE_SQLITE_BLOBS_STAGED(
        IMPORT_QUERIES.out.db,
        PREDICT_STRUCTURES.out.cif.map { _meta, cifs -> cifs },
        ANNOTATE_REPS.out.s4pred_features.map { _meta, dir -> dir }.collect().ifEmpty([]),
        ANNOTATE_REPS.out.tm_features.map { _meta, dir -> dir }.collect().ifEmpty([]),
        PARSE_BIOMES.out.res.map { _meta, files -> files }.collect().ifEmpty([]),
        PARSE_DOMAINS.out.res.map { _meta, files -> files }.collect().ifEmpty([]),
        UPDATE_FAMILIES.out.successful_ids.map { ids -> ids.sort().join('\n') }.collectFile(name: 'successful_ids.txt', newLine: true),
        ch_update_info
    )
    ch_versions = ch_versions.mix( UPDATE_SQLITE_BLOBS_STAGED.out.versions )

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
    ch_multiqc_config        = channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = multiqc_config ?
        channel.fromPath(multiqc_config, checkIfExists: true) :
        channel.empty()
    ch_multiqc_logo          = multiqc_logo ?
        channel.fromPath(multiqc_logo, checkIfExists: true) :
        channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
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
    ch_multiqc_files = ch_multiqc_files.mix(UPDATE_FAMILIES.out.seqkit_stats_mqc.collect { t -> t[1] }.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:
    MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
}
