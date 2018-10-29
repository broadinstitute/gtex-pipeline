task bamsync {

    File source_bam
    File target_bam
    File target_bam_index
    String prefix

    Int memory
    Int disk_space
    Int num_threads
    Int num_preempt

    command {
        set -euo pipefail
        echo $(date +"[%b %d %H:%M:%S] Running bamsync")
        /src/run_bamsync.sh ${source_bam} ${target_bam} ${prefix}
        echo $(date +"[%b %d %H:%M:%S] Running samtools flagstat")
        samtools flagstat ${prefix}.Aligned.sortedByCoord.out.patched.bam > ${prefix}.Aligned.sortedByCoord.out.patched.bam.flagstat.txt
    }

    output {
        File patched_bam_file="${prefix}.Aligned.sortedByCoord.out.patched.bam"
        File patched_bam_index="${prefix}.Aligned.sortedByCoord.out.patched.bam.bai"
        File patched_bam_flagstat="${prefix}.Aligned.sortedByCoord.out.patched.bam.flagstat.txt"
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


workflow bamsync_workflow {
    call bamsync
}
