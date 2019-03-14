task leafcutter_bam_to_junc {

    File bam_file
    String sample_id

    Int memory
    Int disk_space
    Int num_threads
    Int num_preempt

    command {
        set -euo pipefail
        PATH=/opt/leafcutter/scripts:$PATH
        echo $(date +"[%b %d %H:%M:%S] LeafCutter: extracting junctions for sample ${sample_id}")
        # modified /opt/leafcutter/scripts/bam2junc.sh to incorporate WASP filters:
        # sh /opt/leafcutter/scripts/bam2junc.sh ${bam_file} ${sample_id}.leafcutter.junc
        bed_file=${bam_file}.bed
        junc_file=${sample_id}.leafcutter.junc
        samtools view ${bam_file} | grep -v "vW:i:[2-7]" | filter_cs.py | sam2bed.pl --use-RNA-strand - $bed_file
        bed2junc.pl $bed_file $junc_file
        gzip $junc_file
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
        File junc_file="${sample_id}.leafcutter.junc.gz"
    }

    meta {
        author: "Francois Aguet"
    }
}

workflow leafcutter_bam_to_junc_workflow {
    call leafcutter_bam_to_junc
}
