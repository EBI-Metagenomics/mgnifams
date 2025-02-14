/*
    SEQUENCE CLUSTERING
*/

include { EXTRACT_UNANNOTATED_FASTA } from "../../../subworkflows/local/extract_unannotated_fasta"
include { EXECUTE_CLUSTERING        } from "../../../subworkflows/local/execute_clustering"

workflow SETUP_CLUSTERS {
    take:
    input
    fasta_input_mode
    compress_mode

    main:
    if (!fasta_input_mode) {
        mgnifams_input_fa = EXTRACT_UNANNOTATED_FASTA( input, compress_mode ).fasta_ch
    } else {
        mgnifams_input_fa = channel.fromPath(input)
    }
    EXECUTE_CLUSTERING(mgnifams_input_fa)

    emit:
    mgnifams_input_fa = mgnifams_input_fa
    clusters_tsv      = EXECUTE_CLUSTERING.out.clusters_tsv
    num_sequences     = EXECUTE_CLUSTERING.out.num_sequences
}
