task torus {

    File qtl_file
    File annotation_file
    String prefix

    Int memory
    Int disk_space
    Int num_threads
    Int num_preempt

    command {
        set -euo pipefail
        python3 /src/run_torus.py ${qtl_file} ${annotation_file} ${prefix}
    }

    runtime {
        docker: "gcr.io/broad-cga-francois-gtex/torus:latest"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }

    output {
        File torus_output="${prefix}.torus_enrichment.txt"
        File torus_log="${prefix}.torus_enrichment.log"
    }

    meta {
        author: "Francois Aguet"
    }
}


workflow torus_workflow {
    call torus
}
