#!/usr/bin/env nextflow

include { FOLDSEEK_EASYSEARCH as FOLDSEEK_EASYSEARCH_PDB        } from "$launchDir/modules/foldseek/easysearch/main.nf"
include { FOLDSEEK_EASYSEARCH as FOLDSEEK_EASYSEARCH_ALPHAFOLDB } from "$launchDir/modules/foldseek/easysearch/main.nf"
include { FOLDSEEK_EASYSEARCH as FOLDSEEK_EASYSEARCH_ESM        } from "$launchDir/modules/foldseek/easysearch/main.nf"

workflow ANNOTATE_STRUCTURES {
    take:
    pdb_ch
    
    main:
    def foldseek_pdb_path       = params.foldseek_db_path + "/pdb"
    def foldseek_alphafold_path = params.foldseek_db_path + "/alphafold"
    def foldseek_esm_path       = params.foldseek_db_path + "/esm"

    pdb_aln       = Channel.empty()
    alphafold_aln = Channel.empty()
    esm_aln       = Channel.empty()
    
    pdb_db = [ [ id:'pdb' ], file(foldseek_pdb_path) ]
    pdb_aln = FOLDSEEK_EASYSEARCH_PDB(pdb_ch, pdb_db).aln
    if (workflow.profile == "slurm" || workflow.profile == "lsf") {
        alphafold_db = [ [ id:'alphafold' ], file(foldseek_alphafold_path) ]
        alphafold_aln = FOLDSEEK_EASYSEARCH_ALPHAFOLDB(pdb_ch, alphafold_db).aln
        esm_db = [ [ id:'esm' ], file(foldseek_esm_path) ]
        esm_aln = FOLDSEEK_EASYSEARCH_ESM(pdb_ch, esm_db).aln
    }

    emit:
    pdb_aln
    alphafold_aln
    esm_aln
}
