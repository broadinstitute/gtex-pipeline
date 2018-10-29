task samtools_view {

    File bam_file
    File bam_index
    String prefix
    String? options
    String? region
    File? reference_fasta
    File? reference_fasta_index

    Int memory
    Int disk_space
    Int num_threads
    Int num_preempt

    command {
        set -euo pipefail
        echo $(date +"[%b %d %H:%M:%S] Running 'samtools view'.")
        samtools view ${options} ${"-T " + reference_fasta} ${bam_file} ${region} > ${prefix}.bam
        samtools index ${prefix}.bam
        echo $(date +"[%b %d %H:%M:%S] done.")
    }

    output {
        File output_bam_file = "${prefix}.bam"
        File output_bam_index = "${prefix}.bam.bai"
    }

    runtime {
        docker: "gcr.io/broad-cga-francois-gtex/gtex_rnaseq:V9"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }

    meta {
        author: "Francois Aguet"
    }
}


workflow samtools_view_workflow {
    call samtools_view
}
