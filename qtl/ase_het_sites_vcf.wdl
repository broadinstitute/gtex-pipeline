task ase_het_sites_vcf {

    File vcf_file
    String individual_id

    Int memory
    Int disk_space
    Int num_threads
    Int num_preempt

    command {
        set -euo pipefail
        bash /src/get_het_sites_vcf.sh ${vcf_file} ${individual_id}
    }

    output {
        File het_vcf = "${individual_id}.hets.vcf.gz"
        File het_vcf_index = "${individual_id}.hets.vcf.gz.tbi"
    }

    runtime {
        docker: "gcr.io/broad-cga-francois-gtex/gtex_eqtl:V8"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }

    meta {
        author: "Francois Aguet"
    }
}


workflow ase_het_sites_vcf_workflow {
    call ase_het_sites_vcf
}
