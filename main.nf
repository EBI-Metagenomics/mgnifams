#!/usr/bin/env nextflow

include { initiate_proteins                               } from "$baseDir/subworkflows/initiate_proteins/main.nf"
include { execute_clustering                              } from "$baseDir/subworkflows/execute_clustering/main.nf"
include { create_families                                 } from "$baseDir/subworkflows/create_families/main.nf"
include { FASTA_DOMAINANNOTATION                          } from "$baseDir/subworkflows/fasta_domainannotation/main.nf"
include { EXTRACT_UNIQUE_IDS as EXTRACT_UNIQUE_BLASTP_IDS } from "$baseDir/modules/general.nf"
include { EXTRACT_UNIQUE_IDS as EXTRACT_UNIQUE_IPS_IDS    } from "$baseDir/modules/general.nf"
include { FIND_UNANNOTATED_IDS                            } from "$baseDir/modules/general.nf"
// include { EXPORT_INTERPRO_ANNOTATIONS_CSV } from "$baseDir/modules/exporting.nf"
// include { EXPORT_BLASTP_ANNOTATIONS_CSV   } from "$baseDir/modules/exporting.nf"
// include { CONCAT_ANNOTATIONS              } from "$baseDir/modules/exporting.nf"
// include { produce_models                  } from "$baseDir/subworkflows/produce_models/main.nf"
include { annotate_structures                             } from "$baseDir/subworkflows/annotate_structures/main.nf"

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

    blastp_ids = EXTRACT_UNIQUE_BLASTP_IDS(annotations_ch.blastp_tsv
            .map { meta, path -> path }
            .collect())
    interproscan_ids = EXTRACT_UNIQUE_IPS_IDS(annotations_ch.inteproscan_tsv.collect())

    // blastp_csv = EXPORT_BLASTP_ANNOTATIONS_CSV(annotations_ch.blastp_tsv.map { meta, path -> path }).collectFile(name: 'blastp.csv')
    // inteproscan_csv = EXPORT_INTERPRO_ANNOTATIONS_CSV(annotations_ch.inteproscan_tsv.map { meta, path -> path }).collectFile(name: 'interproscan.csv')
    
    annotations_ch = blastp_ids
        .concat(interproscan_ids)
        .collectFile()

    // annotations_ch = CONCAT_ANNOTATIONS(annotations_ch, 'strategy90')
    FIND_UNANNOTATED_IDS(annotations_ch, families.reps_ids)
    // // unknown_models = produce_models(families.families_folder)

    annotate_structures(FIND_UNANNOTATED_IDS.out, families.reps_fasta)
    
}