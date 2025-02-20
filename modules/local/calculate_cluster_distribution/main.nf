process CALCULATE_CLUSTER_DISTRIBUTION {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(meta), path(clusters)

    output:
    tuple val(meta), path("${prefix}_clustering_distribution_mqc.csv"), emit: mqc
    path "versions.yml"                                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    cat <<EOF > ${prefix}_clustering_distribution_mqc.csv
    # id: "cluster_distribution"
    # section_name: "Initial Clustering Sizes Distribution"
    # description: "The seed clusters were calculated via the mmseqs suite. The table below presents cluster sizes and how many times each size was encountered."
    # format: "csv"
    # plot_type: "table"
    Sample Name,Cluster Size,Number of Clusters
    EOF

    awk '{count[\$1]++} END {for (c in count) size[count[c]]++} END {for (s in size) print s "," s "," size[s]}' ${clusters} | sort -n --parallel=${task.cpus} >> ${prefix}_clustering_distribution_mqc.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sed: \$(sed --version 2>&1 | sed -n 1p | sed 's/sed (GNU sed) //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_clustering_distribution_mqc.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}
