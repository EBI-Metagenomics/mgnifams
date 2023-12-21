process PARSE_CIF {
    label "venv"
    publishDir "${params.outdir}", mode: "copy"

    input:
    tuple val(meta), path(pdb_folder)

    output:
    tuple val(meta), path("cif/${meta.id}"), emit: cif

    script:
    """
    mkdir -p cif/${meta.id}
    
    for file in ${pdb_folder}/*; do
        name=\$(basename \$file .pdb)
        python3 ${params.scriptDir}/parse_cif.py \$file cif/${meta.id}/\$name.cif
    done
    """
}
