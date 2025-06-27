#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { RUN_ESMFOLD_CONDA } from '../../run_esmfold_conda/main'

workflow test_esmfold {

    input = [
        [ id:'test' ],
        file("tests/modules/ebi-metagenomics/esmfold/data/one_protein.fasta", checkIfExists: true) // params.test_data['sarscov2']['genome']['proteome_fasta']
    ]

    RUN_ESMFOLD_CONDA ( input, "cpu" )
}

workflow test_esmfold_gpu {

    input = [
        [ id:'test' ],
        file("tests/modules/ebi-metagenomics/esmfold/data/one_protein.fasta", checkIfExists: true)
    ]

    RUN_ESMFOLD_CONDA ( input, "gpu" )
}
