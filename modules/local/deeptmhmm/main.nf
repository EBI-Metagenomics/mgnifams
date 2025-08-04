process DEEPTMHMM {
    tag "$meta.id"
    label 'process_single'
    stageInMode 'copy'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://interpro/deeptmhmm:1.0' :
        'docker.io/interpro/deeptmhmm:1.0' }"

    input:
    tuple val(meta), path(fasta)
    path tmhmm_dir

    output:
    tuple val(meta), path("outdir")

    script:
    """
    # deeptmhmm has a hard coded assumption it is being run within its dir
    cd ${tmhmm_dir}
    python3 predict.py \\
        --fasta ../${fasta} \\
        --output-dir ../outdir
    cd ..
    rm -r ${tmhmm_dir} outdir/embeddings
    chmod -R 777 outdir
    """
}