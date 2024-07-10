process PARSE_CIF {
    label "venv"
    publishDir "${params.outDir}/structures/cif", mode: "copy"

    input:
    tuple val(meta), path(pdb_folder)

    output:
    tuple val(meta), path("*.cif"), emit: cif

    script:
    """    
    for file in ${pdb_folder}/*; do
        name=\$(basename \$file .pdb)
        parse_cif.py \$file \$name.cif
    done
    """
}
