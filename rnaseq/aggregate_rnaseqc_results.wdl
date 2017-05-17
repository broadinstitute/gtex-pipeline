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
        # workaround for broken 'write_lines'
        echo $(date +"[%b %d %H:%M:%S] Writing inputs to file")
        python3 <<CODE
with open("rpkm_gct_paths.tsv", "w") as f:
    for p in '${sep="," rnaseqc_rpkm_gcts}'.split(","):
        f.write(p+'\n')
with open("count_gct_paths.tsv", "w") as f:
    for p in '${sep="," rnaseqc_count_gcts}'.split(","):
        f.write(p+'\n')
with open("exon_count_gct_paths.tsv", "w") as f:
    for p in '${sep="," rnaseqc_exon_count_gcts}'.split(","):
        f.write(p+'\n')
CODE
        echo $(date +"[%b %d %H:%M:%S] Combining RPKM GCTs")
        python3 /src/combine_GCTs.py rpkm_gct_paths.tsv "${prefix}.rnaseqc_rpkm"
        echo $(date +"[%b %d %H:%M:%S] Combining count GCTs")
        python3 /src/combine_GCTs.py count_gct_paths.tsv "${prefix}.rnaseqc_counts"
        echo $(date +"[%b %d %H:%M:%S] Combining exon count GCTs")
        python3 /src/combine_GCTs.py exon_count_gct_paths.tsv "${prefix}.rnaseqc_exon_counts"
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
