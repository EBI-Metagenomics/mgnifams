#!/usr/bin/env nextflow

include { PRODUCE_MODELS } from "../../subworkflows/produce_models/main.nf"

workflow {
    Channel
        .fromPath(params.families_path) 
        .set { families_ch }

    def unannotated_ids = new HashSet()
        file(params.unannotated_ids_path).eachLine { line ->
            line = line.trim()
            if (line.length() > 0) {
                unannotated_ids.add(line)
            }
        }

    families_ch
        .filter { file -> unannotated_ids.contains(file.getBaseName()) }
        .set { unannotated_families_ch }
    
    PRODUCE_MODELS(unannotated_families_ch)
}
