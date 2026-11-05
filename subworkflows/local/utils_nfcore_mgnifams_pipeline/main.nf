//
// Subworkflow with functionality specific to the ebi-metagenomics/mgnifams pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    _monochrome_logs  // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    _input            //  string: Path to input samplesheet

    main:

    ch_versions = channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Create channel from input file provided through params.input
    //
    def samplesheet_rows = samplesheetToList(params.input, "${projectDir}/assets/schema_input.json")

    if (params.mode == "run_mgnifams_pipeline" && !params.hhdb_path) {
        error("run_mgnifams_pipeline: --hhdb_path is required")
    }
    if (params.mode in ["update_mgnifams", "post_update_mgnifams_update_db", "annotate_not_updated", "rebuild_biomes"]) {
        // All collected outputs and the delta DB are per run
        if (samplesheet_rows.size() != 1) {
            error("${params.mode}: the samplesheet must have exactly one row, found ${samplesheet_rows.size()}")
        }
    }
    if (params.mode == "update_mgnifams") {
        if (params.parquet_chunks < 1 || params.hmm_chunk_size < 1) {
            error("update_mgnifams: --parquet_chunks and --hmm_chunk_size must be >= 1")
        }
        if (params.run_alphafold2 && !workflow.stubRun && !params.colabfold_params_path) {
            error("update_mgnifams: --run_alphafold2 needs --colabfold_params_path")
        }
        // Optional reps_fasta / structures (the last two schema columns): re-annotate the not-updated families
        def row = samplesheet_rows[0]
        if (!row[11] != !row[12]) {
            error("update_mgnifams: set both reps_fasta and structures, or neither")
        }
        if (!row[11]) {
            log.warn("update_mgnifams: no reps_fasta/structures, not-updated families keep their previous annotations")
        }
        // The annotation-only merge replaces all fold rows of the families; AlphaFoldDB is searched only under slurm
        else if (!workflow.profile.contains("slurm") && !workflow.profile.contains("test")) {
            error("update_mgnifams: the AlphaFoldDB Foldseek search runs only under a slurm profile; without it the merge would drop the not-updated families' AlphaFoldDB hits")
        }
    }

    // Rows follow the schema_input.json property order
    ch_samplesheet = channel.fromList(samplesheet_rows)

    if (params.mode == "run_mgnifams_pipeline") {
        ch_samplesheet = ch_samplesheet
            .map { sample, protein_input, _results_folder, _existing_db, _schema, _mgnprotein_db_config, _sequences, _clusters, _pfam, _hmms, _families, _reps, _structures, _refined_families, _biome_results ->
                [ sample, file(protein_input, checkIfExists: true) ]
            }
    }
    else if (params.mode == "init_mgnifams_db") {
        ch_samplesheet = ch_samplesheet
            .map { sample, _protein_input, results_folder, _existing_db, schema, _mgnprotein_db_config, _sequences, _clusters, _pfam, _hmms, _families, _reps, _structures, _refined_families, _biome_results ->
                [ sample, schema, results_folder ]
            }
    }
    else if (params.mode == "update_mgnifams_db") {
        ch_samplesheet = ch_samplesheet
            .map { sample, _protein_input, results_folder, existing_db, _schema, mgnprotein_db_config, _sequences, _clusters, _pfam, _hmms, _families, _reps, _structures, _refined_families, _biome_results ->
                [ sample, results_folder, existing_db, mgnprotein_db_config ]
            }
    }
    else if (params.mode == "update_mgnifams") {
        ch_samplesheet = ch_samplesheet
            .map { sample, _protein_input, _results_folder, _existing_db, _schema, mgnprotein_db_config, sequences, clusters, pfam, hmms, families, _reps, _structures, _refined_families, _biome_results ->
                if (!sequences || !clusters || !pfam || !hmms || !families) {
                    error("update_mgnifams: samplesheet row '${sample.id}' needs mgnify_proteins_sequences, mgnify_proteins_clusters, mgnify_proteins_pfam, mgnifams_hmms and mgnifams_families")
                }
                if (mgnprotein_db_config) {
                    log.warn("update_mgnifams: mgnprotein_db_config is deprecated and ignored; biomes come from mgnify_proteins_clusters (cluster_rep_biomes)")
                }
                [ sample, file(sequences), file(clusters), file(pfam), file(hmms), file(families) ]
            }
    }
    else if (params.mode == "annotate_not_updated") {
        // The annotation-only merge replaces all fold rows of the families; AlphaFoldDB is searched only under slurm
        if (!workflow.profile.contains("slurm") && !workflow.profile.contains("test")) {
            error("annotate_not_updated: the AlphaFoldDB Foldseek search runs only under a slurm profile; without it the merge would drop the families' AlphaFoldDB hits")
        }
        ch_samplesheet = ch_samplesheet
            .map { sample, _protein_input, _results_folder, _existing_db, _schema, _mgnprotein_db_config, _sequences, _clusters, _pfam, _hmms, families, reps, structures, _refined_families, _biome_results ->
                if (!families || !reps || !structures) {
                    error("annotate_not_updated: samplesheet row '${sample.id}' needs mgnifams_families, reps_fasta and structures")
                }
                [ sample, file(families), file(reps), file(structures) ]
            }
    }
    else if (params.mode == "rebuild_biomes") {
        // Its pipeline_info and biome_results would overwrite the update's own
        if (file("${params.outdir}/update_families/update_info.csv").exists()) {
            error("rebuild_biomes: --outdir ${params.outdir} is an update_mgnifams outdir; use a separate outdir")
        }
        ch_samplesheet = ch_samplesheet
            .map { sample, _protein_input, _results_folder, _existing_db, _schema, _mgnprotein_db_config, _sequences, clusters, _pfam, _hmms, _families, _reps, _structures, refined_families, _biome_results ->
                if (!refined_families || !clusters) {
                    error("rebuild_biomes: samplesheet row '${sample.id}' needs refined_families and mgnify_proteins_clusters")
                }
                [ sample, file(refined_families), file(clusters) ]
            }
    }
    else if (params.mode == "post_update_mgnifams_update_db") {
        ch_samplesheet = ch_samplesheet
            .map { sample, _protein_input, results_folder, _existing_db, _schema, _mgnprotein_db_config, _sequences, _clusters, _pfam, _hmms, _families, _reps, _structures, _refined_families, biome_results ->
                if (!results_folder) {
                    error("post_update_mgnifams_update_db: samplesheet row '${sample.id}' needs results_folder (an update_mgnifams outdir)")
                }
                // Optional biome_results (a rebuild_biomes outdir's biome_results/) replaces results_folder/biome_results
                [ sample, results_folder, biome_results ?: [] ]
            }
    }

    // update_mgnifams: [ meta, reps fasta or [], structures folder or [] ]; empty for the other modes
    ch_reannotate = params.mode == "update_mgnifams" ?
        channel.fromList(samplesheet_rows).map { row -> [ row[0], row[11] ? file(row[11]) : [], row[12] ? file(row[12]) : [] ] } :
        channel.empty()

    emit:
    samplesheet = ch_samplesheet
    reannotate  = ch_reannotate
    versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications
    multiqc_report  //    value: List of MultiQC report paths

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.getVal(),
            )
        }

        completionSummary(monochrome_logs)
        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Validate channels from input samplesheet
//
//
// Number of families (one rep per family) in a release's families.tsv.gz (data rows, header excluded): hmmsearch -Z for the reps searches
//
def countFamilies(families) {
    return families.withInputStream { stream ->
        new java.util.zip.GZIPInputStream(stream).newReader().readLines().findAll { line -> line.trim() }.size() - 1
    }
}

def validateInputSamplesheet(input) {
    def (metas, fastqs) = input[1..2]

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect{ meta -> meta.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    return [ metas[0], fastqs ]
}
//
// Generate methods description for MultiQC
//
def toolCitationText() {
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def citation_text = [
            "Tools used in the workflow included:",
            "MultiQC (Ewels et al. 2016)",
            "."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def reference_text = [
            "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>"
        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    // Only uncomment below if logic in toolCitationText/toolBibliographyText has been filled!
    // meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    // meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
