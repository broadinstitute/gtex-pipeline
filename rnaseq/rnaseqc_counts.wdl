task rnaseqc_counts {
    File bam_file
    File bam_index
    File genes_gtf
    File genome_fasta
    File genome_fasta_index

    String prefix
    String? gatk_flags

    Int memory
    Int disk_space
    Int num_threads
    Int num_preempt

    command {
        set -euo pipefail
        touch ${bam_index}
        touch ${genome_fasta_index}
        python3 /src/run_rnaseqc.py ${bam_file} ${genes_gtf} ${genome_fasta} ${prefix}\
            --java /usr/lib/jvm/java-1.7.0-openjdk-amd64/bin/java\
            --memory ${memory}\
            --rnaseqc_flags noDoC strictMode\
            ${" --gatk_flags " + gatk_flags}
    }

    output {
        File gene_rpkm = "${prefix}.gene_rpkm.gct.gz"
        File gene_counts = "${prefix}.gene_reads.gct.gz"
        File exon_counts = "${prefix}.exon_reads.gct.gz"
        File count_metrics = "${prefix}.metrics.tsv"
        File count_outputs = "${prefix}.tar.gz"
    }

    runtime {
        docker: "broadinstitute/gtex_rnaseq:V8"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }

    meta {
        author: "Francois Aguet"
    }
}


workflow rnaseqc_counts_workflow {
    call rnaseqc_counts
}
