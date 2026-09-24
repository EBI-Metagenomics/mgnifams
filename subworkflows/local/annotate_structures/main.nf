#!/usr/bin/env nextflow

include { FOLDSEEK_EASYSEARCH as FOLDSEEK_EASYSEARCH_PDB         } from '../../../modules/nf-core/foldseek/easysearch/main'
include { FOLDSEEK_EASYSEARCH as FOLDSEEK_EASYSEARCH_ALPHAFOLDDB } from '../../../modules/nf-core/foldseek/easysearch/main'

workflow ANNOTATE_STRUCTURES {
    take:
    pdb
    foldseek_pdb_db
    foldseek_alphafold_db
    outdir

    main:
    ch_versions                 = channel.empty()
    ch_alphafold_aln            = channel.empty()

    ch_pdb_db = channel.of(file(foldseek_pdb_db, checkIfExists: true)).map { db -> [ [ id:db.name ], db ] }
    FOLDSEEK_EASYSEARCH_PDB( pdb, ch_pdb_db ).aln
    ch_versions = ch_versions.mix(
        FOLDSEEK_EASYSEARCH_PDB.out.versions_foldseek
            .map { process, tool, version -> "\"${process}\":\n    ${tool}: ${version}" }
    )

    if (workflow.profile.contains("slurm") && !workflow.profile.contains("test")) {
        ch_alphafold_db = channel.of(file(foldseek_alphafold_db, checkIfExists: true)).map { db -> [ [ id:db.name ], db ] }
        ch_alphafold_aln = FOLDSEEK_EASYSEARCH_ALPHAFOLDDB( pdb, ch_alphafold_db ).aln
        ch_versions = ch_versions.mix(
            FOLDSEEK_EASYSEARCH_ALPHAFOLDDB.out.versions_foldseek
                .map { process, tool, version -> "\"${process}\":\n    ${tool}: ${version}" }
        )
    }

    ch_foldseek_hits = FOLDSEEK_EASYSEARCH_PDB.out.aln
        .concat(ch_alphafold_aln)
        .map { _meta, file ->
            file
        }
        .collectFile(name: 'all_hits.tsv', keepHeader: true)
    ch_foldseek_hits.collectFile(name: 'all_hits.tsv', storeDir: outdir + "/annotation/structures/foldseek") // // Published copy only: consumers use the work-dir file, so deleting outdir keeps the -resume cache
    ch_foldseek_hits = ch_foldseek_hits
        .map { file ->
            [[id: 'foldseek_hits'], file]
        }

    emit:
    versions      = ch_versions
    foldseek_hits = ch_foldseek_hits
}
