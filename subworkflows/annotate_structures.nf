#!/usr/bin/env nextflow

include { FOLDSEEK_EASYSEARCH as FOLDSEEK_EASYSEARCH_PDB        } from "../modules/nf-core/foldseek/easysearch/main.nf"
include { FOLDSEEK_EASYSEARCH as FOLDSEEK_EASYSEARCH_ALPHAFOLDB } from "../modules/nf-core/foldseek/easysearch/main.nf"
include { FOLDSEEK_EASYSEARCH as FOLDSEEK_EASYSEARCH_ESM        } from "../modules/nf-core/foldseek/easysearch/main.nf"
include { FIND_ANNOTATED_FAMILIES_BY_STRUCTURE                  } from "../modules/local/find_annotated_families_by_structure.nf"

workflow ANNOTATE_STRUCTURES {
    take:
    pdb_ch
    
    main:
    def foldseek_pdb_path       = params.foldseek_db_path + "/pdb"
    def foldseek_alphafold_path = params.foldseek_db_path + "/alphafold"
    def foldseek_esm_path       = params.foldseek_db_path + "/esm"

    alphafold_aln = Channel.empty()
    esm_aln       = Channel.empty()
    
    pdb_db = Channel.of([ [ id:'pdb' ], file(foldseek_pdb_path) ])
    pdb_aln = FOLDSEEK_EASYSEARCH_PDB(pdb_ch, pdb_db).aln
    if (workflow.profile.contains("slurm")) {
        alphafold_db = Channel.of([ [ id:'alphafold' ], file(foldseek_alphafold_path) ])
        alphafold_aln = FOLDSEEK_EASYSEARCH_ALPHAFOLDB(pdb_ch, alphafold_db).aln
        esm_db = Channel.of([ [ id:'esm' ], file(foldseek_esm_path) ])
        esm_aln = FOLDSEEK_EASYSEARCH_ESM(pdb_ch, esm_db).aln
    }

    pdb_aln
        .concat(alphafold_aln)
        .concat(esm_aln)
        .map { meta, file ->
            file
        }
        .collectFile(name: 'foldseek_hits.tsv', storeDir: params.outdir + "/structures/foldseek")
        .map { file ->
            [[id: 'foldseek_hits'], file]
        }
        .set { foldseek_hits }

    FIND_ANNOTATED_FAMILIES_BY_STRUCTURE(foldseek_hits)

    emit:
    foldseek_hits
}
