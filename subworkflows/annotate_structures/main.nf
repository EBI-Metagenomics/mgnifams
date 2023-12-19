#!/usr/bin/env nextflow

include { FOLDSEEK_EASYSEARCH } from "$launchDir/modules/foldseek/easysearch/main.nf" // as FOLDSEEK_EASYSEARCH_PDB

workflow ANNOTATE_STRUCTURES {
    take:
    pdb_ch
    
    main:
    def foldseek_pdb_path = params.foldseek_db_path + "/pdb"
    
    pdb_db = [ [ id:'pdb' ], file(foldseek_pdb_path) ]
    pdb_aln = FOLDSEEK_EASYSEARCH(pdb_ch, pdb_db).aln

    emit:
    pdb_aln
}
