#!/usr/bin/env nextflow

include { initiate_proteins               } from "$baseDir/subworkflows/initiate_proteins/main.nf"
include { execute_clustering              } from "$baseDir/subworkflows/execute_clustering/main.nf"
include { create_families                 } from "$baseDir/subworkflows/create_families/main.nf"
include { FASTA_DOMAINANNOTATION          } from "$baseDir/subworkflows/fasta_domainannotation/main.nf"
include { EXPORT_INTERPRO_ANNOTATIONS_CSV } from "$baseDir/modules/exporting.nf"
include { EXPORT_EGGNOG_ANNOTATIONS_CSV   } from "$baseDir/modules/exporting.nf"
include { EXPORT_BLASTP_ANNOTATIONS_CSV   } from "$baseDir/modules/exporting.nf"
include { CONCAT_ANNOTATIONS              } from "$baseDir/modules/exporting.nf"
include { FIND_UNANNOTATED_IDS            } from "$baseDir/modules/general.nf"
// include { produce_models } from "$baseDir/subworkflows/produce_models/main.nf"

workflow {
    combined_fasta_file = initiate_proteins()
    mmseqs = execute_clustering(combined_fasta_file)
    families = create_families(combined_fasta_file, mmseqs.clu_tsv)
    
    // FASTA_DOMAINANNOTATION -----
    fasta_chunks_ch = families.reps_fasta.splitFasta( by: params.chunk_size, file: true )
    fasta_chunks_ch
        .map { filepath ->
            def parts = filepath.baseName.split('\\.')
            def number = parts[1]
            def id = "mg${number}"
            return [ [id:id], [file(filepath)] ]
        }
        .set { input }

    def uniprot_sprot_fasta_path = params.dataDir + params.uniprot_sprot_fasta_name
    blast_fasta = Channel.of( [ [id:'test'], uniprot_sprot_fasta_path ] )
    annotations_ch = FASTA_DOMAINANNOTATION( input, blast_fasta, "diamond" )
    // ----------------------------

    // blastp_csv = EXPORT_BLASTP_ANNOTATIONS_CSV(annotations_ch.blastp_csv.map { meta, path -> path })
    // inteproscan_csv = EXPORT_INTERPRO_ANNOTATIONS_CSV(annotations_ch.inteproscan_tsv.map { meta, path -> path })
    // annotations_ch = blastp_csv.concat(inteproscan_csv).collect()
    // annotations_ch = CONCAT_ANNOTATIONS(annotations_ch, 'strategy90')
    // FIND_UNANNOTATED_IDS(annotations_ch, families.reps_ids)
    // // unknown_models = produce_models(families.families_folder)
}