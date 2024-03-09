task peer_factors {

    File phenotype_file
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
        Rscript /src/run_PEER.R ${phenotype_file} ${prefix} ${num_peer}
        /src/combine_covariates.py ${prefix}.PEER_covariates.txt ${prefix} ${"--genotype_pcs " + genotype_pcs} ${"--add_covariates " + add_covariates}
        gzip *.PEER_residuals.txt
    }

    runtime {
        docker: "gcr.io/broad-cga-francois-gtex/gtex_eqtl:V10"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }

    output {
        File combined_covariates="${prefix}.combined_covariates.txt"
        File alpha="${prefix}.PEER_alpha.txt"
        File residuals="${prefix}.PEER_residuals.txt.gz"
    }

    meta {
        author: "Francois Aguet"
    }
}

workflow peer_factors_workflow {
    call peer_factors
}
