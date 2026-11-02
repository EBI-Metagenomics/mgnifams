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
        HHSUITE_HHBLITS( HHSUITE_REFORMAT.out.msa, ch_hhdb.first() )
    } else if (hh_mode == "hhsearch") {
        HHSUITE_HHSEARCH( HHSUITE_REFORMAT.out.msa, ch_hhdb.first() )
    } else {
        throw new Exception("Invalid hh_mode value. Should be 'hhblits' or 'hhsearch'.")
    }

    emit:
    hh_mode == "hhblits" ? HHSUITE_HHBLITS.out.hhr : HHSUITE_HHSEARCH.out.hhr
}
