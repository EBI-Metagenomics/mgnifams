#!/usr/bin/env nextflow

include { UNZIP_BZ2     } from "$launchDir/modules/general.nf"
include { REMOVE_HEADER } from "$launchDir/modules/general.nf"

workflow PREPROCESS_INPUT {
    main:
    Channel
        .fromPath(params.mgy90_path)
        .set { mgy90_file_bz2 }

    mgy90_file_with_header = UNZIP_BZ2(mgy90_file_bz2)    
    preprocessed_mgy90_file = REMOVE_HEADER(mgy90_file_with_header)

    emit:
    preprocessed_mgy90_file
}
