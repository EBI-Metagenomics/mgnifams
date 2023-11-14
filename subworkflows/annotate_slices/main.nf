#!/usr/bin/env nextflow

include { INTERPROSCAN } from "$baseDir/modules/interproscan.nf"
include { EGGNOG_MAPPER } from "$baseDir/modules/eggnog.nf"
include { MAKEBLASTDB; BLASTP } from "$baseDir/modules/blast.nf"
include { EXPORT_INTERPRO_ANNOTATIONS_CSV; EXPORT_EGGNOG_ANNOTATIONS_CSV; EXPORT_BLASTP_ANNOTATIONS_CSV; CONCAT_ANNOTATIONS } from "$baseDir/modules/exporting.nf"

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
        .fromPath(uniprot_sprot_fasta_path) 
        .set { uniprot_sprot_fasta_file }

    fasta_chunks_ch = fasta_ch.splitFasta( by: params.chunk_size, file: true )
    interpro_ch = INTERPROSCAN(fasta_chunks_ch).collectFile(name: "ips_annotations.tsv")
    eggnog_ch = EGGNOG_MAPPER(fasta_chunks_ch, ch_eggnong_data_dir.first(), ch_eggnog_diamond_db.first(), ch_eggnog_db.first()).collectFile(name: "eggnog_annotations.tsv")
    library_ch = MAKEBLASTDB(uniprot_sprot_fasta_file)
    blastp_ch = BLASTP(fasta_chunks_ch, library_ch.first()).collectFile(name: "blastp_annotations.tsv")

    interpro_csv = EXPORT_INTERPRO_ANNOTATIONS_CSV(interpro_ch)
    eggnog_csv = EXPORT_EGGNOG_ANNOTATIONS_CSV(eggnog_ch)
    blastp_csv = EXPORT_BLASTP_ANNOTATIONS_CSV(blastp_ch)
    annotations_ch = CONCAT_ANNOTATIONS(interpro_csv.concat(eggnog_csv, blastp_csv).collect(), 'unknown')

    emit:
    annotations_ch
}