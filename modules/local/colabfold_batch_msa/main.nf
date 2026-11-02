process COLABFOLD_BATCH_MSA {
    tag "$meta.id"
    label 'process_medium'
    label 'process_gpu'

    // Adapted from nf-core/proteinfold modules/local/colabfold_batch: predicts from precomputed a3m MSAs
    container "nf-core/proteinfold_colabfold:2.0.0"

    input:
    tuple val(meta), path(a3ms, stageAs: 'input/*') // <family_id>.a3m.gz, representative first
    path colabfold_params, stageAs: 'params/*'      // AlphaFold2 weights (params_model_*.npz); [] only in stub runs
    val num_recycles

    output:
    tuple val(meta), path("*.pdb")        , emit: pdb
    tuple val(meta), path("*_scores.json"), emit: scores
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("Local COLABFOLD_BATCH_MSA module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }
    def args = task.ext.args ?: ''
    """
    set -euo pipefail

    mkdir -p params a3m out
    # A staged weights directory (e.g. alphafold_params_2022-12-06/) is flattened, as ColabFold reads params/ only
    for d in params/*/; do
        [ -d "\$d" ] && ln -sf "\$(realpath "\$d")"/* params/
    done
    compgen -G "params/params_model_*.npz" > /dev/null || { echo "COLABFOLD_BATCH_MSA: no params_model_*.npz in the ColabFold params" >&2; exit 1; }
    touch params/download_finished.txt

    for f in input/*.a3m.gz; do
        gzip -cdf "\$f" > "a3m/\$(basename "\$f" .gz)"
    done

    colabfold_batch \\
        ${args} \\
        --num-recycle ${num_recycles} \\
        --data \$PWD \\
        a3m/ \\
        out/

    for f in a3m/*.a3m; do
        id=\$(basename "\$f" .a3m)
        pdb=( out/"\${id}"_unrelaxed_rank_001_*.pdb )
        scores=( out/"\${id}"_scores_rank_001_*.json )
        [ -f "\${pdb[0]}" ] && [ -f "\${scores[0]}" ] || { echo "COLABFOLD_BATCH_MSA: no rank_001 prediction for \${id}" >&2; exit 1; }
        cp "\${pdb[0]}" "\${id}.pdb"
        cp "\${scores[0]}" "\${id}_scores.json"
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        colabfold_batch: \$(pip list | grep "^colabfold" | awk '{print \$2}' 2>/dev/null || echo "unknown")
    END_VERSIONS
    """

    stub:
    """
    for f in input/*.a3m.gz; do
        id=\$(basename "\$f" .a3m.gz)
        touch "\${id}.pdb" "\${id}_scores.json"
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        colabfold_batch: \$(pip list | grep "^colabfold" | awk '{print \$2}' 2>/dev/null || echo "unknown")
    END_VERSIONS
    """
}
