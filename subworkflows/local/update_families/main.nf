include { EXTRACT_UNANNOTATED_PARQUET_SLICES } from '../../../modules/local/extract_unannotated_parquet_slices/main'
include { CHECK_QUALITY                      } from '../check_quality'
include { HMMER_ESLSFETCHINDEX               } from '../../../modules/nf-core/hmmer/eslsfetchindex/main'
include { SPLIT_HMM_LIB                      } from '../../../modules/local/split_hmm_lib/main'
include { MGNIFAM_UPDATEFAMILIES             } from '../../../modules/nf-core/mgnifam/updatefamilies/main'
include { POOL_UPDATED_FAMILIES              } from '../../../modules/local/pool_updated_families/main'

workflow UPDATE_FAMILIES {
    take:
    ch_samplesheet      // channel: [ meta, sequences_parquet, pfam_parquet, hmm_lib, db_config or [] ]
    n_chunks            // integer: parallel preprocessing tasks
    min_sequence_length // integer
    hmm_chunk_size      // integer: family HMMs per update task
    outdir

    main:
    ch_versions = channel.empty()
    ch_meta     = ch_samplesheet.map { meta, _sequences, _pfam, _hmms, _db_config -> meta }

    // Unannotated slices of the MGnify proteins, one task per group of sequence row groups
    ch_extract_input = ch_samplesheet
        .combine(channel.of(0..<n_chunks))
        .map { meta, sequences, pfam, _hmms, _db_config, chunk_index -> [ meta, sequences, pfam, chunk_index ] }
    EXTRACT_UNANNOTATED_PARQUET_SLICES( ch_extract_input, n_chunks, min_sequence_length )
    ch_versions = ch_versions.mix( EXTRACT_UNANNOTATED_PARQUET_SLICES.out.versions )

    // Zero-padded chunk names: name order is chunk order, so the FASTA is byte-identical across runs
    ch_fasta = EXTRACT_UNANNOTATED_PARQUET_SLICES.out.fa
        .map { _meta, fa -> fa }
        .collectFile(name: 'mgnifams_update.faa', storeDir: "${outdir}/update_families", sort: { fa -> fa.name })
        .combine(ch_meta)
        .map { fa, meta -> [ meta, fa ] }

    CHECK_QUALITY( ch_fasta )
    ch_versions = ch_versions.mix( CHECK_QUALITY.out.versions )

    HMMER_ESLSFETCHINDEX( ch_fasta )
    ch_versions = ch_versions.mix(
        HMMER_ESLSFETCHINDEX.out.versions_hmmer.mix( HMMER_ESLSFETCHINDEX.out.versions_easel )
            .map { process, tool, version -> "\"${process}\":\n    ${tool}: ${version}" }
    )

    SPLIT_HMM_LIB( ch_samplesheet.map { meta, _sequences, _pfam, hmms, _db_config -> [ meta, hmms ] }, hmm_chunk_size )
    ch_versions = ch_versions.mix( SPLIT_HMM_LIB.out.versions )

    // One update task per HMM chunk, all against the same indexed FASTA
    ch_update_input = SPLIT_HMM_LIB.out.hmms
        .flatMap { meta, libs -> [libs].flatten().collect { lib -> [ meta + [chunk: lib.name.tokenize('._')[1]], lib ] } }
        .combine( ch_fasta.join( HMMER_ESLSFETCHINDEX.out.ssi ).map { _meta, fa, ssi -> [ fa, ssi ] } )
    MGNIFAM_UPDATEFAMILIES( ch_update_input )
    ch_versions = ch_versions.mix(
        MGNIFAM_UPDATEFAMILIES.out.versions_mgnifam
            .map { process, tool, version -> "\"${process}\":\n    ${tool}: ${version}" }
    )

    POOL_UPDATED_FAMILIES(
        ch_meta,
        collectedFiles( MGNIFAM_UPDATEFAMILIES.out.tsv ),
        collectedFiles( MGNIFAM_UPDATEFAMILIES.out.csv ),
        collectedFiles( MGNIFAM_UPDATEFAMILIES.out.reps_fasta ),
        collectedFiles( MGNIFAM_UPDATEFAMILIES.out.delta ),
        collectedFiles( MGNIFAM_UPDATEFAMILIES.out.discarded ),
        collectedFiles( MGNIFAM_UPDATEFAMILIES.out.converged )
    )
    ch_versions = ch_versions.mix( POOL_UPDATED_FAMILIES.out.versions )

    // Successful family ids = family_ids.fasta headers; an all-discard update is fatal, after the pool is published
    ch_successful_ids = POOL_UPDATED_FAMILIES.out.family_ids_fasta
        .map { _meta, fa ->
            def ids = [] as Set
            fa.eachLine { line -> if (line.startsWith('>')) ids << line.substring(1) }
            if (!ids) {
                error("update_mgnifams: no successful families; see ${outdir}/update_families/updated_delta.csv")
            }
            ids
        }

    // Full MSAs are per family ([id: family], sto.gz); every successful family must have one
    ch_full_msa = MGNIFAM_UPDATEFAMILIES.out.full_msa
        .flatMap { _meta, msas -> [msas].flatten() }
        .map { msa -> [ [id: msa.name.replace('.sto.gz', '')], msa ] }
        .unique { t -> t[0].id } // ponytail: ids are unique across chunks (SPLIT_HMM_LIB); only the module stub repeats 1_7
    ch_checked_ids = ch_successful_ids
        .map { ids -> [ 'ids', ids ] }
        .join( ch_full_msa.map { t -> t[0].id }.collect().ifEmpty([]).map { ids -> [ 'ids', ids as Set ] } )
        .map { _key, ids, msa_ids ->
            def missing = ids - msa_ids
            if (missing) {
                error("update_mgnifams: successful families without a full MSA: ${missing.sort().take(20).join(', ')}")
            }
            ids
        }
    ch_successful_full_msa = ch_full_msa
        .combine( ch_checked_ids )
        .filter { meta, _msa, ids -> meta.id in ids }
        .map { meta, msa, _ids -> [ meta, msa ] }

    emit:
    family_ids_fasta = POOL_UPDATED_FAMILIES.out.family_ids_fasta.map { _meta, fa -> [ [id: 'reps_fasta'], fa ] }
    successful_ids   = ch_checked_ids                    // value: Set of family ids
    refined_families = POOL_UPDATED_FAMILIES.out.tsv
    metadata         = POOL_UPDATED_FAMILIES.out.metadata.map { _meta, csv -> [ [id: 'metadata'], csv ] }
    delta            = POOL_UPDATED_FAMILIES.out.delta
    full_msa         = ch_successful_full_msa            // [ [id: family], full_msa.sto.gz ]
    seqkit_stats_mqc = CHECK_QUALITY.out.seqkit_stats_mqc
    versions         = ch_versions
}

// All files of an optional per-chunk output, or [] when no chunk produced one
def collectedFiles(ch) {
    return ch.map { _meta, f -> f }.collect().ifEmpty([])
}
