task combine_gcts {

    Array[File] gct_list
    String prefix

    Int memory
    Int disk_space

    command {
        python3 /src/combine_GCTs.py "${sep=' ' gct_list}" ${prefix}
    }

    output {
        File combined_gct="${prefix}.gct.gz"
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


workflow combine_gcts_workflow {
    call combine_gcts
}
