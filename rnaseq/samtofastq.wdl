task samtofastq {

    File input_bam_cram
    String prefix
    File? reference_fasta

    Float memory
    Int java_memory = floor(memory - 0.5)
    Int disk_space
    Int num_threads
    Int num_preempt

    command {
        set -euo pipefail

        # make sure path is absolute
        input_bam_abs=${input_bam_cram}
        if [[ $input_bam_abs != /* ]]; then
            input_bam_abs=$PWD/$input_bam_abs
        fi

        mkdir samtofastq  # workaround for named pipes
        python3 -u /src/run_SamToFastq.py $input_bam_abs -p ${prefix} ${"--reference_fasta " + reference_fasta} --output_dir samtofastq --memory ${java_memory}
        mv samtofastq/${prefix}_*.fastq.gz .
    }

    output {
        File fastq1="${prefix}_1.fastq.gz"
        File fastq2="${prefix}_2.fastq.gz"
    }

    runtime {
        docker: "gcr.io/broad-cga-francois-gtex/gtex_rnaseq:V9"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }

    meta {
        author: "Francois Aguet"
    }
}


workflow samtofastq_workflow {
    call samtofastq
}
