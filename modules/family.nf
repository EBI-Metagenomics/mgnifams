process EXPORT_REPS {
    publishDir 'data/output', mode: 'copy'

    input:
    path clust_tsv

    output:
    path "rep_names.txt"

    script:
    """
    awk -F'\t' '{print \$1}' ${clust_tsv} | sort -u > rep_names.txt
    """
}

process CREATE_FAMILY_FA {    
    input:
    path clust_tsv
    path fasta
    val mgyp

    output:
    path "${mgyp}_family.fa"

    script:
    """
    python3 ${baseDir}/bin/family_rep_into_fasta.py ${clust_tsv} ${fasta} ${mgyp}_family.fa ${mgyp}
    """
}

process KEEP_UNKNOWN {
    publishDir 'data/output/families/unknown', mode: 'copy', saveAs: { filename ->
        def newFilename = filename.replaceAll("_unknown", "")
        "${newFilename}"
    }
    
    input:
    path fasta

    output:
    path "${fasta.baseName}_unknown.fa", optional: true

    script:
    """
    python3 ${baseDir}/bin/keep_unknown_family.py ${fasta} ${fasta.baseName}_unknown.fa 1
    """
}

process KEEP_KNOWN {
    publishDir 'data/output/families/known', mode: 'copy', saveAs: { filename ->
        def newFilename = filename.replaceAll("_known", "")
        "${newFilename}"
    }
    
    input:
    path fasta

    output:
    path "${fasta.baseName}_known.fa", optional: true

    script:
    """
    python3 ${baseDir}/bin/keep_unknown_family.py ${fasta} ${fasta.baseName}_known.fa 0
    """
}