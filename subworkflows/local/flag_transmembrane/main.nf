#!/usr/bin/env nextflow

include { EXTRACT_FIRST_STOCKHOLM_SEQUENCES } from "../../../modules/local/extract_first_stockholm_sequences/main"
include { DEEPTMHMM                         } from "../../../modules/nf-core/deeptmhmm/main"
include { FLAG_TM                           } from "../../../modules/local/flag_tm/main"

workflow FLAG_TRANSMEMBRANE {
    take:
    tm_msa_ch

    main:
    ch_versions = Channel.empty()

    extracted_res = EXTRACT_FIRST_STOCKHOLM_SEQUENCES(tm_msa_ch)
    // TODO ch_versions = ch_versions.mix( EXTRACT_FIRST_STOCKHOLM_SEQUENCES.out.versions )

    prob_ids = extracted_res.prob_ids
    extracted_res.fa
        .map { meta, filepath ->
            [ meta, filepath.splitFasta( by: params.deeptmhmm_chunk_size, file: true ) ]
        }
        .transpose()
        .map { meta, filepath ->
            [ [id: meta.id, chunk: filepath.getBaseName().split('\\.')[-1]], filepath ]
        }
        .set { fa_ch }

    deeptmhmm_res = DEEPTMHMM(fa_ch)
    ch_versions = ch_versions.mix( DEEPTMHMM.out.versions )

    deeptmhmm_res.line3
        .map { meta, file -> file }
        .collectFile(name: 'predicted_topologies.3line', storeDir: params.outdir + "/redundancy/tm")
    gff3_ch   = deeptmhmm_res.gff3
    tm_ids_ch = FLAG_TM(gff3_ch, params.tm_fraction_threshold)
    // TODO ch_versions = ch_versions.mix( FLAG_TM.out.versions )

    tm_ids_ch
        .map { meta, file -> file }
        .collectFile(name: 'tm_ids.txt', storeDir: params.outdir + "/redundancy/tm")
        .map { filepath ->
            return [ [id:"tm_ids"], file(filepath) ]
        }
        .set { tm_ids_ch }

    fa_ch
        .map { meta, filepath ->
            filepath
        }
        .collectFile(name: 'family_reps.fasta')
        .map { filepath ->
            [ [id:"fa_reps"], file(filepath) ]
        }
        .set { fa_ch }

    emit:
    versions  = ch_versions
    tm_ids_ch = tm_ids_ch
    prob_ids  = prob_ids
    fa_ch     = fa_ch
}
