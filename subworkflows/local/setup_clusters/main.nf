include { PIGZ_UNCOMPRESS                    } from '../../../modules/nf-core/pigz/uncompress/main'
include { EXTRACT_UNANNOTATED_PARQUET_SLICES } from '../../../modules/local/extract_unannotated_parquet_slices/main'
include { MASK_MGNIFAM_REGIONS               } from '../../../subworkflows/local/mask_mgnifam_regions'
include { CHECK_QUALITY                      } from '../../../subworkflows/local/check_quality'
include { EXECUTE_CLUSTERING                 } from '../../../subworkflows/local/execute_clustering'
include { CALCULATE_CLUSTER_DISTRIBUTION     } from '../../../modules/local/calculate_cluster_distribution/main'
include { EXTRACT_UNIQUE_CLUSTER_REPS        } from '../../../modules/local/extract_unique_cluster_reps/main'
include { CHUNK_CLUSTERS                     } from '../../../modules/local/chunk_clusters/main'

workflow SETUP_CLUSTERS {
    take:
    input               // channel: [ meta, protein_input FASTA or [], sequences parquet or [], clusters parquet, pfam parquet, mgnifam_regions dir or [] ]
    parquet_chunks      // integer: parallel Pfam-slicing tasks over the sequence parquet row groups
    min_sequence_length
    outdir
    minimum_members
    clusters_chunk_size

    main:
    ch_versions = channel.empty()

    // One plain FASTA for everything downstream (clustering, the SSI index, generation): a pre-sliced FASTA as is
    // (.gz decompressed), or the Pfam-free slices of the sequence parquet, as update_mgnifams makes them
    ch_branch = input.branch { _meta, fasta, _sequences, _clusters, _pfam, _regions ->
        gz     : fasta && fasta.name.endsWith('.gz')
        fasta  : fasta
        parquet: true
    }

    PIGZ_UNCOMPRESS( ch_branch.gz.map { meta, fasta, _sequences, _clusters, _pfam, _regions -> [ meta, fasta ] } )

    ch_slice_input = ch_branch.parquet
        .combine( channel.of(0..<parquet_chunks) )
        .map { meta, _fasta, sequences, clusters, pfam, _regions, chunk_index -> [ meta, sequences, clusters, pfam, chunk_index ] }
    EXTRACT_UNANNOTATED_PARQUET_SLICES( ch_slice_input, parquet_chunks, min_sequence_length )
    ch_versions = ch_versions.mix( EXTRACT_UNANNOTATED_PARQUET_SLICES.out.versions )

    // Zero-padded chunk names: name order is chunk order, so the FASTA is byte-identical across runs.
    // shortcut: one parquet row per run (as update_mgnifams); several rows would be concatenated into one FASTA
    ch_sliced_file = EXTRACT_UNANNOTATED_PARQUET_SLICES.out.fa
        .map { _meta, fa -> fa }
        .collectFile( name: 'mgnifams_input.faa', sort: { fa -> fa.name } )
    ch_sliced_file.collectFile( name: 'mgnifams_input.faa', storeDir: "${outdir}/setup_clusters" ) // Published copy only: consumers use the work-dir file, so deleting outdir keeps the -resume cache
    ch_sliced = ch_sliced_file
        .combine( ch_branch.parquet.map { meta, _fasta, _sequences, _clusters, _pfam, _regions -> meta } )
        .map { fa, meta -> [ meta, fa ] }

    ch_unmasked_fa = ch_branch.fasta
        .map { meta, fasta, _sequences, _clusters, _pfam, _regions -> [ meta, fasta ] }
        .mix( PIGZ_UNCOMPRESS.out.file, ch_sliced )

    // MGnifam regions of earlier releases cut out as well, when the row has a mgnifam_regions folder
    ch_mask = ch_unmasked_fa
        .join( input.map { meta, _fasta, _sequences, _clusters, _pfam, regions -> [ meta, regions ] } )
        .branch { _meta, _fa, regions ->
            mask : regions
            plain: true
        }
    ch_region_files = ch_mask.mask
        .map { meta, _fa, regions ->
            def tsvs = files("${regions}/*.tsv")
            if (!tsvs) {
                error("run_mgnifams_pipeline: mgnifam_regions ${regions} holds no *.tsv region file")
            }
            [ meta, tsvs ]
        }
    // 256 buckets: ~10M regions and ~2M sequences each for release 2.0 (2.6B regions, 508M sequences)
    MASK_MGNIFAM_REGIONS( ch_mask.mask.map { meta, fa, _regions -> [ meta, fa ] }, ch_region_files, 256, min_sequence_length )
    ch_versions = ch_versions.mix( MASK_MGNIFAM_REGIONS.out.versions )
    MASK_MGNIFAM_REGIONS.out.fasta
        .map { _meta, fa -> fa }
        .collectFile( name: 'mgnifams_input.fa', storeDir: "${outdir}/setup_clusters" ) // Published copy only, as above

    ch_mgnifams_input_fa = ch_mask.plain
        .map { meta, fa, _regions -> [ meta, fa ] }
        .mix( MASK_MGNIFAM_REGIONS.out.fasta )

    CHECK_QUALITY( ch_mgnifams_input_fa )
    ch_versions = ch_versions.mix( CHECK_QUALITY.out.versions )

    EXECUTE_CLUSTERING( ch_mgnifams_input_fa )

    CALCULATE_CLUSTER_DISTRIBUTION( EXECUTE_CLUSTERING.out )
    ch_versions = ch_versions.mix( CALCULATE_CLUSTER_DISTRIBUTION.out.versions )

    EXTRACT_UNIQUE_CLUSTER_REPS( EXECUTE_CLUSTERING.out, minimum_members )
    ch_versions = ch_versions.mix( EXTRACT_UNIQUE_CLUSTER_REPS.out.versions )

    ch_cluster_reps_chunks = EXTRACT_UNIQUE_CLUSTER_REPS.out.reps
        .splitText(file:true, by: clusters_chunk_size)
        .map { meta, file ->
            [[id: meta.id, chunk: file.getBaseName(1).split('\\.')[-1]], file]
        }

    CHUNK_CLUSTERS( ch_cluster_reps_chunks, EXECUTE_CLUSTERING.out.first() )
    ch_versions = ch_versions.mix( CHUNK_CLUSTERS.out.versions )

    emit:
    versions          = ch_versions
    mgnifams_input_fa = ch_mgnifams_input_fa
    seqkit_stats_mqc  = CHECK_QUALITY.out.seqkit_stats_mqc
    cluster_distr_mqc = CALCULATE_CLUSTER_DISTRIBUTION.out.mqc
    cluster_chunks    = CHUNK_CLUSTERS.out.tsv
}
