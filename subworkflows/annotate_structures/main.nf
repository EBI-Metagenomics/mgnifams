#!/usr/bin/env nextflow

include { FIND_FASTA_BY_ID    } from "$baseDir/modules/general.nf"
include { ESMFOLD             } from "$baseDir/modules/esmfold/main.nf"
include { FOLDCOMP_COMPRESS   } from "$baseDir/modules/foldcomp/compress/main.nf"
include { FOLDSEEK_EASYSEARCH } from "$baseDir/modules/foldseek/easysearch/main.nf" // as FOLDSEEK_EASYSEARCH_PDB

workflow annotate_structures {
    take:
    ids
    fasta
    
    main:
    FIND_FASTA_BY_ID(ids, fasta)
    pdb_fasta_chunks_ch = FIND_FASTA_BY_ID.out.splitFasta( by: params.pdb_chunk_size, file: true )
    pdb_fasta_chunks_ch
        .map { filepath ->
            def parts = filepath.baseName.split('\\.')
            def number = parts[1]
            def id = "pdb${number}"
            return [ [id:id], [file(filepath)] ]
        }
        .set { input }
    ESMFOLD(input)

    ESMFOLD.out.pdb
        .map { id, files -> 
            def dir = files[0].parent
            return [id, dir]
        }
        .set { pdb_dir }
    FOLDCOMP_COMPRESS(pdb_dir)

    def foldseek_pdb_path = params.foldseek_db_path + "pdb"
    input_db = [ [ id:'pdb' ], [ file(foldseek_pdb_path) ] ]
    FOLDSEEK_EASYSEARCH(ESMFOLD.out.pdb, input_db)

    emit:
    FOLDCOMP_COMPRESS.out.fcz
    FOLDSEEK_EASYSEARCH.out.aln
}