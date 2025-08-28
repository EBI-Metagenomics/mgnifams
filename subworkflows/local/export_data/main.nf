include { EXPORT_MGNIFAMS                               } from '../../../modules/local/export_mgnifams/main'
include { FILTER_EXPORT_DOMTBL as FILTER_EXPORT_PFAM    } from '../../../modules/local/filter_export_domtbl/main'
include { FILTER_EXPORT_DOMTBL as FILTER_EXPORT_FUNFAMS } from '../../../modules/local/filter_export_domtbl/main'
include { EXPORT_PFAMS                                  } from '../../../modules/local/export_pfams/main'
include { EXPORT_FOLDS                                  } from '../../../modules/local/export_folds/main'

workflow EXPORT_DATA {
    take:
    family_metadata
    predicted_scores
    composition
    tm_composition
    pfam_domains
    funfam_domains
    query_hmm_length_threshold
    pfam_hits
    foldseek_hits
    outdir

    main:
    ch_versions = Channel.empty()

    EXPORT_MGNIFAMS( family_metadata, predicted_scores, composition, tm_composition )
    ch_versions = ch_versions.mix( EXPORT_MGNIFAMS.out.versions )

    FILTER_EXPORT_PFAM( pfam_domains, query_hmm_length_threshold )
    ch_versions = ch_versions.mix( FILTER_EXPORT_PFAM.out.versions )

    FILTER_EXPORT_FUNFAMS( funfam_domains, query_hmm_length_threshold )
    ch_versions = ch_versions.mix( FILTER_EXPORT_FUNFAMS.out.versions )

    EXPORT_PFAMS( pfam_hits )
    ch_versions = ch_versions.mix( EXPORT_PFAMS.out.versions )

    EXPORT_PFAMS.out.csv
        .map { meta, file ->
            file
        }
        .collectFile(name: "mgnifam_model_pfams.csv", storeDir: outdir + "/table_data/", keepHeader: true)
        .map { file ->
            [ [id: "mgnifam_pfams"], file ]
        }

    EXPORT_FOLDS( foldseek_hits )
    ch_versions = ch_versions.mix( EXPORT_FOLDS.out.versions )

    emit:
    versions = ch_versions
}
