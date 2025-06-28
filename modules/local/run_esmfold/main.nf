process RUN_ESMFOLD {
    tag "$meta.id"
    label 'process_medium'
    label 'process_gpu'

    container "nf-core/proteinfold_esmfold:1.1.1"

    input:
    tuple val(meta), path(fasta)
    path ('./checkpoints/')
    val numRec

    output:
    tuple val(meta), path ("*.pdb")              , emit: pdb
    tuple val(meta), path("${prefix}_scores.txt"), emit: scores
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("Local RUN_ESMFOLD module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.0.3' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    // KR - note: removed the *.pdb -> tmp.pdb, tmp.pdb  -> esmfold.pdb. Why not just take directly?
    // Only one .pdb per ESMFold run
    """
    awk '{if (\$0 ~ /^>/) {gsub(/^>/, "", \$1); split(\$1, id, " "); gsub("/", "_", id[1]); print ">" id[1]} else {print}}' ${fasta} > parsed.fasta

    esm-fold \
        -i parsed.fasta \
        -o \$PWD \
        -m \$PWD \
        --num-recycles ${numRec} \
        $args > ${prefix}_scores.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        esm-fold: $VERSION
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.0.3' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch "test.pdb"
    touch "${prefix}_scores.txt"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        esm-fold: $VERSION
    END_VERSIONS
    """
}
