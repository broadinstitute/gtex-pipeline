task extract_PIRs {

    File bam_file  # also works with CRAM
    File bam_index
    File vcf
    File vcf_index
    String participant_id
    String chr

    File? reference_fasta
    File? reference_fasta_index

    Int memory
    Int disk_space
    Int num_threads
    Int num_preempt

    command {
        set -euo pipefail
        touch ${vcf_index}
        python3 /src/bam_to_pir.py --vcf ${vcf} --bam ${bam_file} --participant_id ${participant_id} ${"--fasta " + reference_fasta} --chr ${chr} --output_dir .
    }

    runtime {
        docker: "gcr.io/broad-cga-francois-gtex/shapeit2:latest"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }

    output {
        File pir_file = "${participant_id}.${chr}.pir"
    }

    meta {
        author: "Francois Aguet"
    }
}


workflow extract_PIRs_workflow {
    File chr_list_file

    Array[String] chr_list = read_lines(chr_list_file)
    scatter (c in chr_list) {
        call extract_PIRs { input: chr=c }
    }
}
