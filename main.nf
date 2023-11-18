#!/usr/bin/env nextflow

include { INITIATE_PROTEINS                               } from "$baseDir/subworkflows/initiate_proteins/main.nf"
include { EXECUTE_CLUSTERING                              } from "$baseDir/subworkflows/execute_clustering/main.nf"
include { CREATE_FAMILIES                                 } from "$baseDir/subworkflows/create_families/main.nf"
include { FASTA_DOMAINANNOTATION                          } from "$baseDir/subworkflows/fasta_domainannotation/main.nf"
include { EXTRACT_UNIQUE_IDS as EXTRACT_UNIQUE_BLASTP_IDS } from "$baseDir/modules/general.nf"
include { EXTRACT_UNIQUE_IDS as EXTRACT_UNIQUE_IPS_IDS    } from "$baseDir/modules/general.nf"
include { FIND_UNANNOTATED_IDS                            } from "$baseDir/modules/general.nf"
include { ANNOTATE_STRUCTURES                             } from "$baseDir/subworkflows/annotate_structures/main.nf"

workflow {
    combined_fasta_file = INITIATE_PROTEINS()
    mmseqs = EXECUTE_CLUSTERING(combined_fasta_file)
    families = CREATE_FAMILIES(combined_fasta_file, mmseqs.clu_tsv)
    
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

    def uniprot_sprot_fasta_path = params.dataDir + "/" + params.uniprot_sprot_fasta_name
    blast_fasta = Channel.of( [ [id:'test'], uniprot_sprot_fasta_path ] )
    annotations_ch = FASTA_DOMAINANNOTATION( input, blast_fasta, "diamond" )
    // ----------------------------

    blastp_ids = EXTRACT_UNIQUE_BLASTP_IDS(annotations_ch.blastp_tsv
            .map { meta, path -> path }
            .collect())
    interproscan_ids = EXTRACT_UNIQUE_IPS_IDS(annotations_ch.inteproscan_tsv.collect())
    annotations_ch = blastp_ids
        .concat(interproscan_ids)
        .collectFile()
    FIND_UNANNOTATED_IDS(annotations_ch, families.reps_ids)
    ANNOTATE_STRUCTURES(FIND_UNANNOTATED_IDS.out, families.reps_fasta)
}
