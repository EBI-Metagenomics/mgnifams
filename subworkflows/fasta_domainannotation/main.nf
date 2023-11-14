include { BLAST_MAKEBLASTDB } from "$baseDir/modules/blast/makeblastdb/main"
include { BLAST_BLASTP      } from "$baseDir/modules/blast/blastp/main"
include { DIAMOND_MAKEDB    } from "$baseDir/modules/diamond/makedb/main"
include { DIAMOND_BLASTP    } from "$baseDir/modules/diamond/blastp/main"
include { INTERPROSCAN      } from "$baseDir/modules/interproscan/main"

workflow FASTA_DOMAINANNOTATION {

    take:
    ch_fasta        // channel: [ val(meta), path(fasta) ]
    val_blast_fasta // value: /path/to/reference/fasta for blast
    val_blast_mode  // value: blast or diamond

    main:

    ch_versions = Channel.empty()

    if (val_blast_mode == "blast") {
        BLAST_MAKEBLASTDB ( val_blast_fasta.map { meta, db -> db } )
        ch_versions = ch_versions.mix(BLAST_MAKEBLASTDB.out.versions)
        BLAST_BLASTP ( ch_fasta, BLAST_MAKEBLASTDB.out.db.first(), 'tsv' )
        ch_versions = ch_versions.mix(BLAST_BLASTP.out.versions)
        blastp_tsv = BLAST_BLASTP.out.tsv
    } else if (val_blast_mode == "diamond") {
        DIAMOND_MAKEDB ( val_blast_fasta )
        ch_versions = ch_versions.mix(DIAMOND_MAKEDB.out.versions)
        blast_columns = '' // defaults: qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore
        DIAMOND_BLASTP ( ch_fasta, DIAMOND_MAKEDB.out.db.map { meta, dmnd -> dmnd }.first(), 'txt', blast_columns )
        ch_versions = ch_versions.mix(DIAMOND_BLASTP.out.versions)
        blastp_tsv = DIAMOND_BLASTP.out.txt
    } else {
        throw new Exception("Invalid mode value '$val_blast_mode'. Should be 'blast' or 'diamond'.")
    }

    INTERPROSCAN ( ch_fasta.map { meta, fasta -> fasta }, params.interproscan_db )

    emit:
    blastp_tsv      = blastp_tsv           // channel: [ val(meta), [ tsv ] ]
    inteproscan_tsv = INTERPROSCAN.out.tsv // channel: [ val(meta), [ tsv ] ]
    versions        = ch_versions          // channel: [ versions.yml ]
}
