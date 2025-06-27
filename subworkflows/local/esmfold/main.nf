include { RUN_ESMFOLD } from '../../../modules/local/run_esmfold'

workflow ESMFOLD {

    take:
    ch_samplesheet    // channel: samplesheet read in from --input
    ch_versions       // channel: [ path(versions.yml) ]
    ch_esmfold_params // directory: /path/to/esmfold/params/
    ch_num_recycles   // int: Number of recycles for esmfold

    main:
    ch_multiqc_files = Channel.empty()

    RUN_ESMFOLD(
        ch_samplesheet,
        ch_esmfold_params,
        ch_num_recycles
    )
    ch_versions = ch_versions.mix(RUN_ESMFOLD.out.versions)

    RUN_ESMFOLD.out.pdb
        .map { meta, file ->
            meta = meta.clone();
            meta.model = "esmfold";
            [ meta, file ]
        }
        .set{ch_pdb_final}

    RUN_ESMFOLD
        .out
        .multiqc
        .map { it[1] }
        .toSortedList()
        .map { [ [ "model": "esmfold"], it.flatten() ] }
        .set { ch_multiqc_report  }

    emit:
    pdb            = ch_pdb_final      // channel: [ id, /path/to/*.pdb ]
    multiqc_report = ch_multiqc_report // channel: /path/to/multiqc_report.html
    versions       = ch_versions       // channel: [ path(versions.yml) ]
}
