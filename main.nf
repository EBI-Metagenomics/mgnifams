#!/usr/bin/env nextflow

include { initiate_proteins      } from "$baseDir/subworkflows/initiate_proteins/main.nf"
include { execute_clustering     } from "$baseDir/subworkflows/execute_clustering/main.nf"
include { create_families        } from "$baseDir/subworkflows/create_families/main.nf"
include { FASTA_DOMAINANNOTATION } from "$baseDir/subworkflows/fasta_domainannotation/main.nf"

include { produce_models         } from "$baseDir/subworkflows/produce_models/main.nf"
include { FIND_UNANNOTATED_IDS   } from "$baseDir/modules/general.nf"

workflow {
    combined_fasta_file = initiate_proteins()
    mmseqs = execute_clustering(combined_fasta_file)
    families = create_families(combined_fasta_file, mmseqs.clu_tsv)

    // FASTA_DOMAINANNOTATION -----
    input = Channel.of( [ [id:'mgnifams_annotate'], families.reps_fasta ] )
    def uniprot_sprot_fasta_path = params.dataDir + params.uniprot_sprot_fasta_name
    blast_fasta = Channel.value( uniprot_sprot_fasta_path )
    eggnog_db = Channel.value( file(params.eggnog_db) ) 
    eggnog_data_dir = eggnog_db.parent
    eggnog_diamond_db = Channel.fromPath(params.eggnong_diamond_db) 
    // annotations_ch = FASTA_DOMAINANNOTATION( input, blast_fasta, eggnog_db, eggnog_data_dir, eggnog_diamond_db)
    // ----------------------------

    // unknown_models = produce_models(families.families_folder)
    // FIND_UNANNOTATED_IDS(annotations_ch, families.reps_ids)
}