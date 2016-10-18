task rsem_reference {

    File reference_fasta
    File annotation_gtf
    String prefix

    Int memory
    Int disk_space
    Int num_threads

    command {
        mkdir ${prefix} && cd ${prefix}
        rsem-prepare-reference ${reference_fasta} rsem_reference --gtf ${annotation_gtf} --num-threads ${num_threads}
        cd .. && tar -cvzf ${prefix}.tar.gz ${prefix}
    }

    output {
        File rsem_reference = "${prefix}.tar.gz"
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


workflow rsem_reference_workflow {
    call rsem_reference
}
