task rnaseqc_metrics {
    File bam_file
    File bam_index
    File genes_gtf
    File genome_fasta
    File genome_fasta_index

    String prefix
    String notes

    Int memory
    Int disk_space
    Int num_preempt

    command {
        touch ${bam_index}
        touch ${genome_fasta_index}
        /usr/lib/jvm/java-1.7.0-openjdk-amd64/bin/java -Xmx${memory}g -jar /opt/RNA-SeQC_1.1.9/RNA-SeQC.jar -n 1000 \
        -s ${prefix},${bam_file},${notes} -t ${genes_gtf} -r ${genome_fasta} -o .
        mv metrics.tsv ${prefix}.metrics.tsv
        mkdir report
        mv report.html report/
        mv meanCoverage* report/
        mv gapLengthHist* report/
        tar -cvzf ${prefix}_report.tar.gz report/*
        tar -cvzf ${prefix}.tar.gz ${prefix}/*
    }

    output {
        File rnaseqc_metrics = "${prefix}.metrics.tsv"
        File rnaseqc_report = "${prefix}_report.tar.gz"
        File rnaseqc_metrics_outputs = "${prefix}.tar.gz"
    }

    runtime {
        docker: "broadinstitute/gtex_rnaseq"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        preemptible: "${num_preempt}"
    }

    meta {
        author: "Francois Aguet"
    }
}


workflow rnaseqc_metrics_workflow {
    call rnaseqc_metrics
}
