include { HHSUITE_REFORMAT } from '../../../modules/nf-core/hhsuite/reformat/main'
include { HHSUITE_HHBLITS  } from '../../../modules/nf-core/hhsuite/hhblits/main'
include { HHSUITE_HHSEARCH } from '../../../modules/nf-core/hhsuite/hhsearch/main'

workflow ANNOTATE_MODELS {
    take:
    seed_msa
    hh_mode
    hhdb_path

    main:
    HHSUITE_REFORMAT( seed_msa, "fas", "a3m" )

    ch_hhdb = channel.of([ [ id: 'pfam_hh_db' ], file(hhdb_path, checkIfExists: true) ])
    if (hh_mode == "hhblits") {
        ch_hhr = HHSUITE_HHBLITS( HHSUITE_REFORMAT.out.msa, ch_hhdb.first() ).hhr
    } else if (hh_mode == "hhsearch") {
        ch_hhr = HHSUITE_HHSEARCH( HHSUITE_REFORMAT.out.msa, ch_hhdb.first() ).hhr
    } else {
        throw new Exception("Invalid hh_mode value. Should be 'hhblits' or 'hhsearch'.")
    }

    emit:
    pfam_hits = ch_hhr
}
