#!/usr/bin/env nextflow

include { EXTRACT_FIRST_STOCKHOLM_SEQUENCES } from "$launchDir/modules/general.nf"
// include { FIND_FASTA_BY_ID  } from "$launchDir/modules/general.nf"
// include { GET_FIRST_SEQUENCE_FROM_STOCKHOLM } from "$launchDir/modules/general.nf"
include { ESMFOLD           } from "$launchDir/modules/esmfold/main.nf"
// include { FOLDCOMP_COMPRESS } from "$launchDir/modules/foldcomp/compress/main.nf"

workflow PREDICT_STRUCTURES {
    take:
    msa
    unannotated_ids
    mode
    
    main:
    family_reps_fasta = EXTRACT_FIRST_STOCKHOLM_SEQUENCES(msa, unannotated_ids, mode)
    family_reps_fasta.view()

    // pdb_fasta_chunks_ch = ids.splitFasta( by: params.pdb_chunk_size, file: true )
    // pdb_fasta_chunks_ch
    //     .map { filepath ->
    //         def parts = filepath.baseName.split('\\.')
    //         def number = parts[1]
    //         def id = "pdb${number}"
    //         return [ [id:id], [file(filepath)] ]
    //     }
    //     .set { input }
    // ESMFOLD(input)

    // ESMFOLD.out.pdb
    //     .map { id, files -> 
    //         def dir = files[0].parent
    //         return [id, dir]
    //     }
    //     .set { pdb_dir }
    // FOLDCOMP_COMPRESS(pdb_dir)

    // def foldseek_pdb_path = params.foldseek_db_path + "pdb"
    // input_db = [ [ id:'pdb' ], [ file(foldseek_pdb_path) ] ]
    // FOLDSEEK_EASYSEARCH(ESMFOLD.out.pdb, input_db)

    // emit:
    // FOLDCOMP_COMPRESS.out.fcz
    // FOLDSEEK_EASYSEARCH.out.aln
}
