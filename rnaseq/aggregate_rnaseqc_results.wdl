task aggregate_rnaseqc_results {

    Array[File] rnaseqc_rpkm_gcts
    Array[File] rnaseqc_count_gcts
    Array[File] rnaseqc_exon_count_gcts
    String prefix

    Int memory
    Int disk_space
    Int num_threads
    Int num_preempt

    command {
        echo $(date +"[%b %d %H:%M:%S] Combining RPKM GCTs")
        python3 /src/combine_GCTs.py ${write_lines(rnaseqc_rpkm_gcts)} "${prefix}.rnaseqc_rpkm"
        echo $(date +"[%b %d %H:%M:%S] Combining count GCTs")
        python3 /src/combine_GCTs.py ${write_lines(rnaseqc_count_gcts)} "${prefix}.rnaseqc_counts"
        echo $(date +"[%b %d %H:%M:%S] Combining exon count GCTs")
        python3 /src/combine_GCTs.py ${write_lines(rnaseqc_exon_count_gcts)} "${prefix}.rnaseqc_exon_counts"
    }

    output {
        File rnaseqc_rpkm_gct="${prefix}.rnaseqc_rpkm.gct.gz"
        File rnaseqc_count_gct="${prefix}.rnaseqc_counts.gct.gz"
        File rnaseqc_exon_count_gct="${prefix}.rnaseqc_exon_counts.gct.gz"
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


workflow aggregate_rnaseqc_results_workflow {
    call aggregate_rnaseqc_results
}
