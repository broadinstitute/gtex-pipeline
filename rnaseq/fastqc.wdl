task fastqc {

    File fastq1
    File? fastq2

    Float memory
    Int disk_space
    Int num_threads
    Int num_preempt

    String fastq1_name = sub(sub(basename(fastq1), "\\.fastq.gz$", ""), "\\.fq.gz$", "" )
    String fastq2_name = sub(sub(basename(fastq2), "\\.fastq.gz$", ""), "\\.fq.gz$", "" )

    command <<<
        set -euo pipefail
        fastqc ${fastq1} ${fastq2} \
            --threads ${num_threads} \
            --outdir .
        unzip -p ${fastq1_name}_fastqc.zip ${fastq1_name}_fastqc/fastqc_data.txt | gzip > ${fastq1_name}.fastqc_data.txt.gz
        unzip -p ${fastq2_name}_fastqc.zip ${fastq2_name}_fastqc/fastqc_data.txt | gzip > ${fastq2_name}.fastqc_data.txt.gz
    >>>

    output {
        File fastq1_fastqc_html = "${fastq1_name}_fastqc.html"
        File fastq1_fastqc_zip =  "${fastq1_name}_fastqc.zip"
        File fastq1_fastqc_data = "${fastq1_name}.fastqc_data.txt.gz"
        File fastq2_fastqc_html = "${fastq2_name}_fastqc.html"
        File fastq2_fastqc_zip =  "${fastq2_name}_fastqc.zip"
        File fastq2_fastqc_data = "${fastq2_name}.fastqc_data.txt.gz"
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

workflow fastqc_workflow {
    call fastqc
}
