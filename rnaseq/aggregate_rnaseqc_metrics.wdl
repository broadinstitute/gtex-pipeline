task aggregate_rnaseqc_metrics {

    Array[File] rnaseqc_sample_metrics
    Array[String] tissue_site_detail
    Array[String] sample_ids
    
    String prefix

    Int memory
    Int disk_space
    Int num_threads
    Int num_preempt

    command {
        # workaround for broken 'write_lines'
        echo $(date +"[%b %d %H:%M:%S] Writing inputs to file")
        python3 <<CODE
with open("metrics_paths.tsv", "w") as f:
    for i,p in zip('${sep="," sample_ids}'.split(","), '${sep="," rnaseqc_sample_metrics}'.split(",")):
        f.write(i+'\t'+p+'\n')
with open("tissue_site_detail.tsv", "w") as f:
    for i,p in zip('${sep="," sample_ids}'.split(","), '${sep="," tissue_site_detail}'.split(",")):
        f.write(i+'\t'+p+'\n')
CODE
        echo $(date +"[%b %d %H:%M:%S] Combining metrics TSVs")
        python3 /src/aggregate_rnaseqc_metrics.py metrics_paths.tsv ${prefix} --annotation_headers TissueSite --annotation_tsvs tissue_site_detail.tsv
    }

    output {
        File rnaseqc_metrics="${prefix}.metrics.tsv"
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


workflow aggregate_rnaseqc_metrics_workflow {
    call aggregate_rnaseqc_metrics
}
