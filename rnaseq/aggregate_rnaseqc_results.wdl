task aggregate_rnaseqc_results {

    Array[File] rnaseqc_count_gcts
    Array[File] rnaseqc_rpkm_gcts
    String prefix

    Int memory
    Int disk_space

    command {
        python3 /src/combine_GCTs.py ${sep=' ' rnaseqc_count_gcts} "${prefix}_rnaseqc_counts"
        python3 /src/combine_GCTs.py ${sep=' ' rnaseqc_rpkm_gcts} "${prefix}_rnaseqc_rpkm"
    }

    output {
        File combined_count_gct="${prefix}_rnaseqc_counts.gct.gz"
        File combined_rpkm_gct="${prefix}_rnaseqc_rpkm.gct.gz"
    }

    runtime {
        docker: "broadinstitute/gtex_rnaseq"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
    }

    meta {
        author: "Francois Aguet"
    }
}


workflow aggregate_rnaseqc_results_workflow {
    call aggregate_rnaseqc_results
}
