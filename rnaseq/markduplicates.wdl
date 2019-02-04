task markduplicates {

    File input_bam
    String prefix
    Int? max_records_in_ram

    Int memory
    Int disk_space
    Int num_threads
    Int num_preempt

    String output_bam = sub(basename(input_bam), "\\.bam$", ".md.bam")

    command {
        set -euo pipefail
        python3 -u /src/run_MarkDuplicates.py ${input_bam} ${prefix} --memory ${memory} ${"--max_records_in_ram " + max_records_in_ram}
        samtools index ${output_bam}
    }

    output {
        File bam_file = "${output_bam}"
        File bam_index = "${output_bam}.bai"
        File metrics = "${prefix}.marked_dup_metrics.txt"
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


workflow markduplicates_workflow {
    call markduplicates
}
