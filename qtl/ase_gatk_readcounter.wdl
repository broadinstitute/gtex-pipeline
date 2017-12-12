task ase_gatk_readcounter {

    File gatk_jar
    File genome_fasta
    File genome_fasta_index
    File genome_fasta_dict
    File het_vcf
    File het_vcf_index
    File bam_file
    File bam_index
    String prefix

    Int memory
    Int disk_space
    Int num_threads
    Int num_preempt

    command {
        set -euo pipefail
        python3 /src/run_GATK_ASEReadCounter.py ${gatk_jar} ${genome_fasta} ${het_vcf} ${bam_file} ${prefix}
    }

    output {
        File ase_read_counts = "${prefix}.readcounts.txt.gz"
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


workflow ase_gatk_readcounter_workflow {
    call ase_gatk_readcounter
}
