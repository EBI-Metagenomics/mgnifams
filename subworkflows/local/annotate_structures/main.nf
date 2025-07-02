#!/usr/bin/env nextflow

include { FOLDSEEK_EASYSEARCH as FOLDSEEK_EASYSEARCH_PDB        } from '../../../modules/nf-core/foldseek/easysearch/main'
include { FOLDSEEK_EASYSEARCH as FOLDSEEK_EASYSEARCH_ALPHAFOLDB } from '../../../modules/nf-core/foldseek/easysearch/main'
include { FIND_ANNOTATED_FAMILIES_BY_STRUCTURE                  } from '../../../modules/local/find_annotated_families_by_structure/main'

workflow ANNOTATE_STRUCTURES {
    take:
    pdb
    foldseek_db_path
    outdir
    
    main:
    ch_versions                 = Channel.empty()
    ch_alphafold_aln            = Channel.empty()
    def foldseek_pdb_path       = foldseek_db_path + "/pdb"
    def foldseek_alphafold_path = foldseek_db_path + "/alphafold"

    ch_pdb_db = Channel.of([ [ id:'pdb' ], file(foldseek_pdb_path, checkIfExists: true) ])
    FOLDSEEK_EASYSEARCH_PDB( pdb, ch_pdb_db ).aln
    ch_versions = ch_versions.mix( FOLDSEEK_EASYSEARCH_PDB.out.versions )

    if (workflow.profile.contains("slurm") && !workflow.profile.contains("test")) {
        ch_alphafold_db = Channel.of([ [ id:'alphafold' ], file(foldseek_alphafold_path, checkIfExists: true) ])
        ch_alphafold_aln = FOLDSEEK_EASYSEARCH_ALPHAFOLDB( pdb, ch_alphafold_db ).aln
        ch_versions = ch_versions.mix( FOLDSEEK_EASYSEARCH_ALPHAFOLDB.out.versions )
    }

    ch_foldseek_hits = FOLDSEEK_EASYSEARCH_PDB.out.aln
        .concat(ch_alphafold_aln)
        .map { meta, file ->
            file
        }
        .collectFile(name: 'all_hits.tsv', storeDir: outdir + "/annotation/structures/foldseek")
        .map { file ->
            [[id: 'foldseek_hits'], file]
        }

    FIND_ANNOTATED_FAMILIES_BY_STRUCTURE( ch_foldseek_hits )
    ch_versions = ch_versions.mix( FIND_ANNOTATED_FAMILIES_BY_STRUCTURE.out.versions )

    emit:
    versions      = ch_versions
    foldseek_hits = ch_foldseek_hits
}
