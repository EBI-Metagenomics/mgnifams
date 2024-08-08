#!/usr/bin/env nextflow

include { EXTRACT_FIRST_STOCKHOLM_SEQUENCES } from "../modules/local/extract_first_stockholm_sequences.nf"
include { DEEPTMHMM                         } from "../modules/nf-core/deeptmhmm/main.nf"
include { FLAG_TM                           } from "../modules/local/flag_tm.nf"

workflow FLAG_TRANSMEMBRANE {
    take:
    tm_msa_ch

    main:
    extracted_res   = EXTRACT_FIRST_STOCKHOLM_SEQUENCES(tm_msa_ch)
    prob_ids        = extracted_res.prob_ids
    fa_ch           = extracted_res.fa
    fa_ch           = fa_ch.map { meta, filepath -> filepath }
    fasta_chunks_ch = fa_ch.splitFasta( by: params.deeptmhmm_chunk_size, file: true )
    fasta_chunks_ch
        .map { filepath ->
            def parts = filepath.baseName.split('\\.')
            def number = parts[1]
            def id = "pdb${number}"
            return [ [id:id], [file(filepath)] ]
        }
        .set { fa_ch }

    deeptmhmm_res = DEEPTMHMM(fa_ch)
    deeptmhmm_res.line3
        .map { it[1] }
        .collectFile(name: 'predicted_topologies.3line', storeDir: params.outdir + "/redundancy/tm")
    gff3_ch   = deeptmhmm_res.gff3
    tm_ids_ch = FLAG_TM(gff3_ch, params.tm_fraction_threshold)
    tm_ids_ch
        .map { it[1] }
        .collectFile(name: 'tm_ids.txt', storeDir: params.outdir + "/redundancy/tm")
        .map { filepath ->
            return [ [id:"tm_ids"], file(filepath) ]
        }
        .set { tm_ids_ch }
    fa_ch
        .map { it[1][0] }
        .collectFile(name: 'family_reps.fasta')
        .map { filepath ->
            return [ [id:"fa_reps"], file(filepath) ]
        }
        .set { fa_ch }

    emit:
    tm_ids_ch = tm_ids_ch
    prob_ids  = prob_ids
    fa_ch     = fa_ch
}
