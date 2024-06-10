process BUILD_HH_DB {
    publishDir "${params.outDir}/hh", mode: "copy"

    input:
    tuple val(meta), path(a3m)

    output:
    tuple val(meta), path("${meta.id}_hh_db"), emit: hh_db

    script:
    """
    ffindex_build -s ${a3m}_a3m.ff{data,index} ${a3m}
    mpirun -np ${task.cpus} ffindex_apply_mpi ${a3m}_a3m.ff{data,index} -i ${a3m}_hhm.ffindex -d ${a3m}_hhm.ffdata -- hhmake -i stdin -o stdout -v 0
    mpirun -np ${task.cpus} cstranslate_mpi -f -x 0.3 -c 4 -I a3m -i ${a3m}_a3m -o ${a3m}_cs219
    sort -k3 -n -r ${a3m}_cs219.ffindex | cut -f1 > sorting.dat
    ffindex_order sorting.dat ${a3m}_hhm.ff{data,index} ${a3m}_hhm_ordered.ff{data,index}
    mv ${a3m}_hhm_ordered.ffindex ${a3m}_hhm.ffindex
    mv ${a3m}_hhm_ordered.ffdata ${a3m}_hhm.ffdata
    ffindex_order sorting.dat ${a3m}_a3m.ff{data,index} ${a3m}_a3m_ordered.ff{data,index}
    mv ${a3m}_a3m_ordered.ffindex ${a3m}_a3m.ffindex
    mv ${a3m}_a3m_ordered.ffdata ${a3m}_a3m.ffdata

    mkdir -p ${meta.id}_hh_db
    mv ${a3m}_* ${meta.id}_hh_db/
    """
}
