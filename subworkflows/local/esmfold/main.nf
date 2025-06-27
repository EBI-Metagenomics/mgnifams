include { RUN_ESMFOLD } from '../../../modules/local/run_esmfold'

workflow ESMFOLD {

    take:
    fasta          // channel: samplesheet read in from --input
    esmfold_params // directory: /path/to/esmfold/params/
    num_recycles   // int: Number of recycles for esmfold

    main:
    ch_versions = Channel.empty()

    RUN_ESMFOLD( fasta, esmfold_params, num_recycles )
    ch_versions = ch_versions.mix( RUN_ESMFOLD.out.versions )

    ch_pdb = RUN_ESMFOLD.out.pdb
        .map { meta, file ->
            meta = meta.clone();
            meta.model = "esmfold";
            [ meta, file ]
        }

    ch_multiqc_report = RUN_ESMFOLD.out.multiqc
        .map { meta, file -> file }
        .toSortedList()
        .map { [ [ "model": "esmfold"], it.flatten() ] }

    emit:
    versions       = ch_versions
    pdb            = ch_pdb            // channel: [ id, /path/to/*.pdb ]
    multiqc_report = ch_multiqc_report // channel: /path/to/multiqc_report.html
}
