task rnaseqc_counts {
    File bam_file
    File bam_index
    File genes_gtf
    File genome_fasta
    File genome_fasta_index

    String prefix
    String notes
    String? gatk_flags

    Int memory
    Int disk_space
    Int num_threads
    Int num_preempt

    command {
        set -euo pipefail
        touch ${bam_index}
        touch ${genome_fasta_index}

        echo $(date +"[%b %d %H:%M:%S] Running RNA-SeQC")
        /usr/lib/jvm/java-1.7.0-openjdk-amd64/bin/java -Xmx${memory}g -jar /opt/RNA-SeQC_1.1.9/RNA-SeQC.jar -n 1000 \
        -s ${prefix},${bam_file},${notes} -t ${genes_gtf} -r ${genome_fasta} -o . -noDoC -strictMode ${" -gatkFlags " + gatk_flags}

        # remove tmp files
        rm genes.rpkm.gct
        rm ${prefix}/${prefix}.metrics.tmp.txt
        rm ${prefix}/${prefix}.metrics.txt
        
        mv ${prefix}/${prefix}.transcripts.rpkm.gct ${prefix}.gene_rpkm.gct
        python3 /src/convert_rnaseqc_counts.py ${prefix}.gene_rpkm.gct ${prefix}/${prefix}.exon_intron_report.txt --exon_report ${prefix}/${prefix}.exon_report.txt ${genes_gtf} -o .
        mv metrics.tsv ${prefix}.metrics.tsv
        gzip ${prefix}.gene_rpkm.gct
        tar -cvzf ${prefix}.tar.gz ${prefix}/*
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
