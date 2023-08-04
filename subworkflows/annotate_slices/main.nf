#!/usr/bin/env nextflow

include { INTERPROSCAN } from "$baseDir/modules/interproscan.nf"
include { EGGNOG_MAPPER } from "$baseDir/modules/eggnog.nf"
include { HMMBUILD; HMMSCAN } from "$baseDir/modules/hmm.nf"
include { EXPORT_INTERPRO_ANNOTATIONS_CSV; EXPORT_EGGNOG_ANNOTATIONS_CSV; EXPORT_UNIPROT_ANNOTATIONS_CSV; CONCAT_ANNOTATIONS } from "$baseDir/modules/exporting.nf"

workflow annotate_slices {
    take:
    fasta_ch
    
    main:
    def uniprot_sprot_fasta_path = params.dataDir + params.uniprot_sprot_fasta_name
    
    Channel
        .fromPath(params.eggnong_data_dir) 
        .set { ch_eggnong_data_dir }
    Channel
        .fromPath(params.eggnong_diamond_db) 
        .set { ch_eggnog_diamond_db }
    Channel
        .fromPath(params.eggnog_db) 
        .set { ch_eggnog_db }
    Channel
        fasta_ch
        .splitFasta(by: 1)
        .view()
        //.set { ch_fasta_files }
    Channel
        .fromPath(uniprot_sprot_fasta_path) 
        .set { uniprot_sprot_fasta_file }

    interpro_ch = INTERPROSCAN(fasta_ch)
    eggnog_ch = EGGNOG_MAPPER(fasta_ch, ch_eggnong_data_dir, ch_eggnog_diamond_db, ch_eggnog_db)
    // build_ch = HMMBUILD(ch_fasta_files)
    interpro_csv = EXPORT_INTERPRO_ANNOTATIONS_CSV(interpro_ch)
    eggnog_csv = EXPORT_EGGNOG_ANNOTATIONS_CSV(eggnog_ch)
    // tblout_ch = HMMSCAN(build_ch, uniprot_sprot_fasta_file.first()).tblout_ch
    uniprot_csv = EXPORT_UNIPROT_ANNOTATIONS_CSV(tblout_ch)
    // CONCAT_ANNOTATIONS(interpro_csv.concat(eggnog_csv, uniprot_csv.collect()).collect(), 'unknown')
}