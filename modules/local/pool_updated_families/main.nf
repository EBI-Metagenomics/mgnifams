process POOL_UPDATED_FAMILIES {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    // Every MGNIFAM_UPDATEFAMILIES output is optional, so each input may be []
    input:
    val meta
    path tsv      , stageAs: 'tsv/*'
    path metadata , stageAs: 'metadata/*'
    path reps     , stageAs: 'reps/*'
    path delta    , stageAs: 'delta/*'
    path discarded, stageAs: 'discarded/*'
    path converged, stageAs: 'converged/*'

    output:
    tuple val(meta), path("refined_families.tsv")  , emit: tsv
    tuple val(meta), path("family_metadata.csv")   , emit: metadata
    tuple val(meta), path("family_reps.fasta")     , emit: reps_fasta
    tuple val(meta), path("family_ids.fasta")      , emit: family_ids_fasta
    tuple val(meta), path("updated_delta.csv")     , emit: delta
    tuple val(meta), path("updated_discarded.csv") , emit: discarded
    tuple val(meta), path("converged_families.txt"), emit: converged
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    set -euo pipefail

    mkdir -p tsv metadata reps delta discarded converged
    # Chunk names are zero-padded, so name order is chunk order
    files() { find -L "\$1" -type f -print0 | sort -z; }

    files tsv | xargs -0 -r cat > refined_families.tsv
    # No header: export_mgnifams.py reads it headerless (same as the run_mgnifams_pipeline family_metadata.csv)
    files metadata | xargs -0 -r -n 1 tail -n +2 > family_metadata.csv
    files reps | xargs -0 -r zcat -f > family_reps.fasta
    sed -E 's/^>[^\\t]*\\t/>/' family_reps.fasta > family_ids.fasta
    files converged | xargs -0 -r cat > converged_families.txt
    {
        echo "family_id,model_length_before,model_length_after,round1_recruits,full_msa_size,retention,rounds_run,converged,outcome"
        files delta | xargs -0 -r -n 1 tail -n +2
    } > updated_delta.csv
    {
        echo "representative,reason,value"
        files discarded | xargs -0 -r -n 1 tail -n +2
    } > updated_discarded.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sed: \$(sed --version 2>&1 | sed -n 1p | sed 's/sed (GNU sed) //')
    END_VERSIONS
    """

    stub:
    // One synthetic successful family, named like the MGNIFAM_UPDATEFAMILIES stub full MSA (1_7)
    """
    printf '1_7\\t1/1-30\\n' > refined_families.tsv
    echo '1_7,1,"1",1-30,30,ACDEFGHIKLMNPQRSTVWYACDEFGHIKL,acdefghiklmnpqrstvwyacdefghikl,False' > family_metadata.csv
    printf '>1/1-30\\t1_7\\nACDEFGHIKLMNPQRSTVWYACDEFGHIKL\\n' > family_reps.fasta
    printf '>1_7\\nACDEFGHIKLMNPQRSTVWYACDEFGHIKL\\n' > family_ids.fasta
    touch converged_families.txt
    printf 'family_id,model_length_before,model_length_after,round1_recruits,full_msa_size,retention,rounds_run,converged,outcome\\n1_7,30,30,1,1,1.0,1,False,successful\\n' > updated_delta.csv
    echo "representative,reason,value" > updated_discarded.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sed: \$(sed --version 2>&1 | sed -n 1p | sed 's/sed (GNU sed) //')
    END_VERSIONS
    """
}
