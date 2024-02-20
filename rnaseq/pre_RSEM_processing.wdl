task remove_IDS_reads {

    File transcriptome_bam
    String prefix

    Int memory
    Int disk_space
    Int num_threads
    Int num_preempt

    command {
        set -euo pipefail
        echo $(date +"[%b %d %H:%M:%S] Running remove_IDS_reads")
        /src/run_remove_IDS_reads.sh ${transcriptome_bam} ${prefix}
    }

    output {
        File transcriptome_noIDS_bam = "${prefix}.Aligned.toTranscriptome_noIDS.out.bam"
    }

    runtime {
        docker: "gcr.io/broad-cga-francois-gtex/gtex_rnaseq:V10"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }

    meta {
        author: "Francois Aguet"
    }
}


workflow pre_RSEM_processing_workflow {
    call remove_IDS_reads
}
