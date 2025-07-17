include { EXPORT_MGNIFAMS } from '../../../modules/local/export_mgnifams/main'
include { EXPORT_FUNFAMS  } from '../../../modules/local/export_funfams/main'
include { EXPORT_PFAMS    } from '../../../modules/local/export_pfams/main'
include { EXPORT_FOLDS    } from '../../../modules/local/export_folds/main'

workflow EXPORT_DATA {
    take:
    family_metadata
    predicted_scores
    composition
    funfam_domains
    pfam_hits
    foldseek_hits
    outdir

    main:
    ch_versions = Channel.empty()

    EXPORT_MGNIFAMS( family_metadata, predicted_scores, composition )
    ch_versions = ch_versions.mix( EXPORT_MGNIFAMS.out.versions )

    EXPORT_FUNFAMS( funfam_domains )
    ch_versions = ch_versions.mix( EXPORT_FUNFAMS.out.versions )

    EXPORT_PFAMS( pfam_hits )
    ch_versions = ch_versions.mix( EXPORT_PFAMS.out.versions )

    EXPORT_PFAMS.out.csv
        .map { meta, file ->
            file
        }
        .collectFile(name: "mgnifam_pfams.csv", storeDir: outdir + "/table_data/", keepHeader: true)
        .map { file ->
            [ [id: "mgnifam_pfams"], file ]
        }

    EXPORT_FOLDS( foldseek_hits )
    ch_versions = ch_versions.mix( EXPORT_FOLDS.out.versions )

    emit:
    versions = ch_versions
}
