#!/usr/bin/env nextflow

include { EXPORT_REPS; CREATE_FAMILY_FA; EXPORT_REPS_FA } from "$baseDir/modules/family.nf"
include { EXPORT_FAMILIES_CSV as EXPORT_KNOWN_FAMILIES_CSV; EXPORT_FAMILIES_CSV as EXPORT_UNKNOWN_FAMILIES_CSV} from "$baseDir/modules/exporting.nf"

workflow create_families {
    take: 
    clust_tsv
    fastaFile

    main:
    reps_file = EXPORT_REPS(clust_tsv)
    reps_ch = reps_file.splitText().map { it.trim() } // removing new line chars at end of mgyps
    families = CREATE_FAMILY_FA(clust_tsv.first(), fastaFile.first(), reps_ch.take(params.debug_top_n_reps)) // TODO, remove take, take all
    
    known_ch = families.known_ch
    known_ch
        .map { file -> file.getBaseName().tokenize('_')[0] }
        .collectFile(name: 'known_ids.txt', newLine: true) 
        .set { known_ids_file }  
    EXPORT_KNOWN_FAMILIES_CSV(known_ids_file, 'known')
    
    unknown_ch = families.unknown_ch
    unknown_reps_fasta = EXPORT_REPS_FA(unknown_ch.collect())
    unknown_ch
        .map { file -> file.getBaseName().tokenize('_')[0] }
        .collectFile(name: 'unknown_ids.txt', newLine: true) 
        .set { unknown_ids_file }  
    EXPORT_UNKNOWN_FAMILIES_CSV(unknown_ids_file, 'unknown')

    emit:
    known_ch
    unknown_ch
    unknown_reps_fasta
}