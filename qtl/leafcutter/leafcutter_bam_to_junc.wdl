task leafcutter_bam_to_junc {

    File bam_file
    String sample_id
    Int? strand_specificity = 0

    Int memory
    Int disk_space
    Int num_threads
    Int num_preempt

    command {
        set -euo pipefail
        echo $(date +"[%b %d %H:%M:%S] Extracting junctions for sample ${sample_id}")
        # select uniquely mapped reads that pass WASP filters
        filtered_bam=${bam_file}.filtered.bam
        samtools view -h -q 255 ${bam_file} | grep -v "vW:i:[2-7]" | samtools view -b > $filtered_bam
        samtools index $filtered_bam
        regtools junctions extract -a 8 -m 50 -M 500000 ${"-s " + strand_specificity} $filtered_bam | gzip -c > ${sample_id}.regtools_junc.txt.gz
        echo $(date +"[%b %d %H:%M:%S] Done")
    }

    runtime {
        docker: "gcr.io/broad-cga-francois-gtex/leafcutter:latest"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }

    output {
        File junc_file="${sample_id}.regtools_junc.txt.gz"
    }

    meta {
        author: "Francois Aguet"
    }
}

workflow leafcutter_bam_to_junc_workflow {
    call leafcutter_bam_to_junc
}
