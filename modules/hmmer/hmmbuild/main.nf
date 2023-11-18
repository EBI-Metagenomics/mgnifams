process HMMER_HMMBUILD {
    publishDir "${params.outdir}/hmm", mode: "copy"
    tag "$meta.id"
    label 'process_low'

    // conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmmer:3.4--hdbdd923_0':
        'biocontainers/hmmer:3.4--hdbdd923_0' }"

    input:
    tuple val(meta), path(alignment)
    path mxfile

    output:
    tuple val(meta), path("*.hmm.gz"), emit: hmm
    // path "*.hmmbuild.txt",             emit: hmmbuildout
    path "versions.yml",               emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}"
    def mxfileopt = mxfile ? "--mxfile ${mxfile}" : ""

    """
    hmmbuild \\
        $args \\
        --cpu $task.cpus \\
        -n ${prefix}  \\
        -o ${prefix}.hmmbuild.txt \\
        ${mxfileopt} \\
        ${prefix}.hmm \\
        $alignment

    gzip ${prefix}.hmm

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hmmer: \$(echo \$(hmmbuild -h | grep HMMER | sed 's/# HMMER //' | sed 's/ .*//' 2>&1))
    END_VERSIONS
    """
}
