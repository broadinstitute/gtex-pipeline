task eqtl_peer_factors {

    File expression_file
    String prefix
    Int num_peer

    File? genotype_pcs
    File? add_covariates

    Int memory
    Int disk_space
    Int num_threads
    Int num_preempt

    command {
        set -euo pipefail
        Rscript /src/run_PEER.R ${expression_file} ${prefix} ${num_peer}
        /src/combine_covariates.py ${prefix}.PEER_covariates.txt ${prefix} ${"--genotype_pcs " + genotype_pcs} ${"--add_covariates " + add_covariates}
    }

    runtime {
        docker: "gcr.io/broad-cga-francois-gtex/gtex_eqtl:V8"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }

    output {
        File combined_covariates="${prefix}.combined_covariates.txt"
        File alpha="${prefix}.PEER_alpha.txt"
    }

    meta {
        author: "Francois Aguet"
    }
}

workflow eqtl_peer_factors_workflow {
    call eqtl_peer_factors
}
