/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ANNOTATE_NOT_UPDATED: re-annotate the not-updated families of a finished update
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    The families.tsv.gz of an update_mgnifams run names them (status not_updated); their
    representatives and structures come from the release that last computed them. Only the
    annotations rerun, with the current databases; --mode post_update_mgnifams_update_db on
    <outdir>/annotate_not_updated builds the annotation-only delta DB.
----------------------------------------------------------------------------------------
*/

include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { countFamilies          } from '../subworkflows/local/utils_nfcore_mgnifams_pipeline'
include { ANNOTATE_FAMILY_SUBSET } from '../subworkflows/local/annotate_family_subset'

workflow ANNOTATE_NOT_UPDATED {

    take:
    ch_samplesheet // channel: [ meta, families.tsv.gz, reps fasta, structures folder ]
    skip_deeptmhmm
    deeptmhmm_path
    pfam_path
    funfams_path
    s4pred_chunk_size
    foldseek_pdb_db
    foldseek_alphafold_db
    query_hmm_length_threshold
    mgnifams_release
    mgnify_proteins_release
    outdir

    main:
    ch_families = ch_samplesheet.map { meta, families, _reps, _structures -> [ meta, families ] }

    ANNOTATE_FAMILY_SUBSET(
        ch_families,
        ch_samplesheet.map { meta, _families, reps, _structures -> [ meta, reps ] },
        ch_samplesheet.map { meta, _families, _reps, structures -> [ meta, structures ] },
        ch_families.map { _meta, families -> countFamilies(families) }.first(),
        skip_deeptmhmm, deeptmhmm_path, pfam_path, funfams_path, s4pred_chunk_size,
        foldseek_pdb_db, foldseek_alphafold_db, query_hmm_length_threshold,
        mgnifams_release, mgnify_proteins_release, outdir
    )

    ANNOTATE_FAMILY_SUBSET.out.n_selected
        .combine( ch_families )
        .subscribe { n, _meta, families ->
            if (n == 0) {
                error("annotate_not_updated: no not_updated families in ${families}")
            }
        }

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

    softwareVersionsToYAML( ANNOTATE_FAMILY_SUBSET.out.versions )
        .mix( topic_versions_string )
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name: 'nf_core_pipeline_software_mqc_versions.yml',
            sort: true,
            newLine: true
        )
}
