include { SPLIT_INTO_BUCKETS as SPLIT_FASTA_BUCKETS   } from '../../../modules/local/mask_mgnifam_regions/main'
include { SPLIT_INTO_BUCKETS as SPLIT_REGIONS_BUCKETS } from '../../../modules/local/mask_mgnifam_regions/main'
include { MASK_BUCKET                                 } from '../../../modules/local/mask_mgnifam_regions/main'

//
// MGnifam regions (every family<TAB>MGYP/start-end row; a bare MGYP = the whole protein) cut out of a protein FASTA:
// both split into hash buckets, each bucket masked in parallel, joined in bucket order (deterministic output).
//
workflow MASK_MGNIFAM_REGIONS {
    take:
    ch_fasta            // channel: [ meta, plain FASTA ] (one row: buckets are joined by number only)
    ch_regions          // channel: [ meta, [ region TSVs ] ]
    buckets             // integer: hash buckets, int(protein id) % buckets
    min_sequence_length // integer: shortest piece kept

    main:
    ch_versions = channel.empty()

    SPLIT_FASTA_BUCKETS( ch_fasta, 'fasta', buckets )
    ch_versions = ch_versions.mix( SPLIT_FASTA_BUCKETS.out.versions )

    SPLIT_REGIONS_BUCKETS( ch_regions, 'regions', buckets )
    ch_versions = ch_versions.mix( SPLIT_REGIONS_BUCKETS.out.versions )

    ch_fasta_buckets = SPLIT_FASTA_BUCKETS.out.buckets
        .flatMap { meta, files -> [files].flatten().collect { f -> [ f.baseName, meta, f ] } }
    ch_region_buckets = SPLIT_REGIONS_BUCKETS.out.buckets
        .flatMap { _meta, files -> [files].flatten().collect { f -> [ f.baseName, f ] } }

    MASK_BUCKET(
        ch_fasta_buckets
            .join( ch_region_buckets, failOnMismatch: true, failOnDuplicate: true )
            .map { bucket, meta, fa, tsv -> [ meta, bucket, fa, tsv ] },
        min_sequence_length
    )
    ch_versions = ch_versions.mix( MASK_BUCKET.out.versions )

    ch_masked = MASK_BUCKET.out.fa
        .map { _meta, fa -> fa }
        .collectFile( name: 'mgnifams_input.fa', sort: { fa -> fa.name } )
        .combine( ch_fasta.map { meta, _fa -> meta } )
        .map { fa, meta -> [ meta, fa ] }

    emit:
    fasta    = ch_masked   // channel: [ meta, masked FASTA ]
    versions = ch_versions
}
