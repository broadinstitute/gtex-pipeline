task rsem_reference {

    File reference_fasta
    File annotation_gtf
    String prefix

    Int memory
    Int disk_space
    Int num_threads
    Int num_preempt

    command {
        mkdir ${prefix} && cd ${prefix}
        rsem-prepare-reference ${reference_fasta} rsem_reference --gtf ${annotation_gtf} --num-threads ${num_threads}
        cd .. && tar -cvzf ${prefix}.tar.gz ${prefix}
    }

    output {
        File rsem_reference = "${prefix}.tar.gz"
    }

    runtime {
        docker: "gcr.io/broad-cga-francois-gtex/gtex_rnaseq:V10"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }

    meta {
        author: "Francois Aguet"
    }
}


workflow rsem_reference_workflow {
    call rsem_reference
}
