/*
    ANNOTATE_FAMILY_SUBSET: re-annotate one status of families (not_updated) of a release from their previous
    representatives and structures: s4pred, DeepTMHMM, Pfam/FunFam hmmsearch (-Z = the release's family count) and
    Foldseek, exported as DB tables. Everything is published under <outdir>/annotate_not_updated/ (conf/modules.config),
    with update_families/update_info.csv (annotation_only,true) for post_update_mgnifams_update_db.
    An empty selection runs nothing.
*/

include { SELECT_FAMILY_SUBSET } from '../../../modules/local/select_family_subset/main'
include { ANNOTATE_REPS        } from '../../../subworkflows/local/annotate_reps'
include { ANNOTATE_STRUCTURES  } from '../../../subworkflows/local/annotate_structures'
include { EXPORT_DATA          } from '../../../subworkflows/local/export_data'

workflow ANNOTATE_FAMILY_SUBSET {
    take:
    families                   // channel: [ meta, families.tsv.gz ]
    reps                       // channel: [ meta, reps fasta ]
    structures                 // channel: [ meta, folder of <family_id>.pdb ]
    z                          // value channel: Integer, the release's family count (hmmsearch -Z)
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
    ch_versions = channel.empty()
    def subset_outdir = "${outdir}/annotate_not_updated"
    // Same conditions as ANNOTATE_REPS (DeepTMHMM) and ANNOTATE_STRUCTURES (AlphaFoldDB)
    def tm_computed   = !skip_deeptmhmm && !workflow.profile.contains("conda")
    def afdb_computed = workflow.profile.contains("slurm") && !workflow.profile.contains("test")

    SELECT_FAMILY_SUBSET( families.join( reps ).join( structures ) )
    ch_versions = ch_versions.mix( SELECT_FAMILY_SUBSET.out.versions )

    ch_fasta = SELECT_FAMILY_SUBSET.out.fasta.filter { _meta, fa -> fa.size() > 0 }
    ch_ids   = ch_fasta.map { _meta, fa -> fa.readLines().findAll { line -> line.startsWith('>') }.collect { line -> line.substring(1) } }

    ANNOTATE_REPS(
        ch_fasta.combine( z ).map { _meta, fa, n -> [ [ id: 'reps_fasta', z: n ], fa ] },
        skip_deeptmhmm, deeptmhmm_path, pfam_path, funfams_path, s4pred_chunk_size
    )
    ch_versions = ch_versions.mix( ANNOTATE_REPS.out.versions )

    ANNOTATE_STRUCTURES( SELECT_FAMILY_SUBSET.out.pdb, foldseek_pdb_db, foldseek_alphafold_db, subset_outdir )
    ch_versions = ch_versions.mix( ANNOTATE_STRUCTURES.out.versions )

    // mgnifam.csv of the selected ids with the recomputed s4pred / TM columns only; the merge reads nothing else
    ch_metadata = ch_ids
        .map { ids -> ( [ 'family_id' ] + ids ).join('\n') }
        .collectFile( name: 'selected_families.csv', newLine: true )
        .map { csv -> [ [ id: 'metadata' ], csv ] }

    EXPORT_DATA( ch_metadata, channel.value([ [ id: 'scores' ], [] ]), \
        ANNOTATE_REPS.out.composition, ANNOTATE_REPS.out.tm_composition, channel.value([ [ id: 'seed_sizes' ], [] ]), \
        ANNOTATE_REPS.out.pfam_domains, ANNOTATE_REPS.out.funfam_domains, query_hmm_length_threshold, \
        channel.empty(), ANNOTATE_STRUCTURES.out.foldseek_hits, subset_outdir )
    ch_versions = ch_versions.mix( EXPORT_DATA.out.versions )

    ch_fasta
        .map { _meta, fa -> fa }
        .collectFile( name: 'family_ids.fasta', storeDir: "${subset_outdir}/update_families" )

    // What this run computed, for the annotation-only delta DB (post_update_mgnifams_update_db)
    ch_fasta
        .combine( z )
        .map { _meta, _fa, n ->
            [
                'key,value',
                "annotation_only,true",
                "tm_computed,${tm_computed}",
                "biome_computed,false",
                "afdb_computed,${afdb_computed}",
                "hmmsearch_z,${n}",
                "pfam_lib,${file(pfam_path).name}",
                "funfams_lib,${file(funfams_path).name}",
                "mgnifams_release,${mgnifams_release}",
                "mgnify_proteins_release,${mgnify_proteins_release}",
                "pipeline_version,${workflow.manifest.version}"
            ].join('\n')
        }
        .collectFile( name: 'update_info.csv', storeDir: "${subset_outdir}/update_families", newLine: true )

    emit:
    n_selected = SELECT_FAMILY_SUBSET.out.fasta.map { _meta, fa -> fa.readLines().count { line -> line.startsWith('>') } }
    versions   = ch_versions
}
