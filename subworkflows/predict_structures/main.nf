#!/usr/bin/env nextflow

include { EXTRACT_FIRST_STOCKHOLM_SEQUENCES_FROM_FOLDER } from "../../modules/family/main.nf"
include { ESMFOLD                                       } from "../../modules/esmfold/main.nf"
include { EXTRACT_LONG_FA                               } from "../../modules/esmfold/extract_long_fa.nf"
include { ESMFOLD_CPU                                   } from "../../modules/esmfold/cpu.nf"
include { EXTRACT_ESMFOLD_SCORES                        } from "../../modules/esmfold/extract_esmfold_scores.nf"
include { PARSE_CIF                                     } from "../../modules/esmfold/parse_cif.nf"

workflow PREDICT_STRUCTURES {
    take:
    msa_sto_ch
    
    main:
    family_reps_fa_ch = EXTRACT_FIRST_STOCKHOLM_SEQUENCES_FROM_FOLDER(msa_sto_ch)
    fasta_chunks_ch = family_reps_fa_ch.splitFasta( by: params.pdb_chunk_size, file: true )
    fasta_chunks_ch
        .map { filepath ->
            def parts = filepath.baseName.split('\\.')
            def number = parts[1]
            def id = "pdb${number}"
            return [ [id:id], [file(filepath)] ]
        }
        .set { fa_ch }
    esmfold_result = ESMFOLD(fa_ch, params.compute_mode)

    // Long sequences that cannot be run on GPU
    fa_ch.map { id, filepath ->
            return filepath
        }
        .set { fa_paths_ch }
    esmfold_result.scores.map { id, filepath ->
            return filepath
        }
        .set { score_paths_ch }

    long_reps_fa_ch = EXTRACT_LONG_FA(fa_paths_ch.collect(), score_paths_ch.collect())
    fasta_long_chunks_ch = long_reps_fa_ch.splitFasta( by: params.pdb_chunk_size, file: true )
    fasta_long_chunks_ch
        .map { filepath ->
            def parts = filepath.baseName.split('\\.')
            def number = parts[1]
            def id = "pdb${number}"
            return [ [id:id], [file(filepath)] ]
        }
        .set { fa_long_ch }
    esmfold_long_result = ESMFOLD_CPU(fa_long_ch)
    // End long sequences

    scores_ch = EXTRACT_ESMFOLD_SCORES(esmfold_result.scores.concat(esmfold_long_result.scores)).csv.map { meta, filepath -> filepath }
    scores_ch = scores_ch.collectFile(name: "pdb_scores.csv", storeDir: params.outdir + "/structures")
    pdb_ch = esmfold_result.pdb.concat(esmfold_long_result.pdb)
    cif_ch = PARSE_CIF(pdb_ch)
    cif_ch
        .map { id, filepath ->
            return filepath }
        .collect()
        .set { cif_ch }

    emit:
    scores_ch
    pdb_ch
    cif_ch
}
