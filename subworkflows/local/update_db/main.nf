include { INSERT_SQLITE } from "../../../modules/local/insert_sqlite/main.nf"
// include { APPEND_BLOBS_PARALLEL } from "../../../modules/local/append_blobs_parallel.nf"

workflow UPDATE_DB {
    take:
    samplesheet
    
    main:

    ch_insertion = samplesheet
        .map { meta, pipeline_results, db ->
            [ meta, file("${pipeline_results.toUriString()}/post-processing/tables/*", checkIfExists: true), db ]
        }
    
    INSERT_SQLITE( ch_insertion )

    // TODO
    //     APPEND_BLOBS_PARALLEL(db, cif_ch, seed_msa_ch, msa_ch, hmm_ch, rf_ch, \
    //         biome_results, domain_results)

}
