version 1.0

# CRAM to FASTQ Conversion Workflow
# 
# This workflow converts CRAM files back to FASTQ format through a three-step process:
# 1. Convert CRAM to aligned BAM using samtools
# 2. Revert aligned BAM to unmapped BAM (uBAM) using GATK RevertSam
# 3. Extract paired-end and unpaired FASTQ files from uBAM using GATK SamToFastq
#
# The workflow outputs compressed FASTQ files (.fastq.gz) for:
# - Read 1 (forward reads)
# - Read 2 (reverse reads) 
# - Unpaired reads (singletons)
#
# This is useful for re-processing sequencing data through different analysis pipelines
# or when you need FASTQ files from archived CRAM data.

workflow cram_to_fastq_workflow {
  meta {
    author: "Shadi Zaheri"
    description: "Convert CRAM files to FASTQ format via BAM and uBAM intermediate steps"
  }
  input {
    File cram_file
    File reference_fasta
    File reference_fasta_index
    File reference_dict
    String sample_id
    String cram_to_bam_memory = "4G"
    Int cram_to_bam_cpu = 2
    String cram_to_bam_disk = "local-disk 20 HDD"
    Int cram_to_bam_preemptible = 1

    String bam_to_ubam_memory = "4G"
    Int bam_to_ubam_cpu = 2
    String bam_to_ubam_disk = "local-disk 20 HDD"
    Int bam_to_ubam_preemptible = 1

    String ubam_to_fastq_memory = "4G"
    Int ubam_to_fastq_cpu = 2
    String ubam_to_fastq_disk = "local-disk 20 HDD"
    Int ubam_to_fastq_preemptible = 1
  }

  call CramToBam {
    input:
      cram = cram_file,
      reference = reference_fasta,
      memory = cram_to_bam_memory,
      cpu = cram_to_bam_cpu,
      disk = cram_to_bam_disk,
      preemptible = cram_to_bam_preemptible
  }

  call BamToUbam {
    input:
      bam = CramToBam.bam,
      reference = reference_fasta,
      memory = bam_to_ubam_memory,
      cpu = bam_to_ubam_cpu,
      disk = bam_to_ubam_disk,
      preemptible = bam_to_ubam_preemptible
  }

  call UbamToFastq {
    input:
      ubam = BamToUbam.ubam,
      memory = ubam_to_fastq_memory,
      cpu = ubam_to_fastq_cpu,
      disk = ubam_to_fastq_disk,
      preemptible = ubam_to_fastq_preemptible,
      prefix = sample_id
  }

  output {
    # Only compressed FASTQ outputs to match task outputs
    File read1_fastq_gz = UbamToFastq.read1_fastq_gz
    File read2_fastq_gz = UbamToFastq.read2_fastq_gz
    File unpaired_fastq_gz = UbamToFastq.unpaired_fastq_gz
  }
}

task CramToBam {
  input {
    File cram
    File reference
    String memory
    Int cpu
    String disk
    Int preemptible
  }
  command {
    samtools view -b -T ~{reference} ~{cram} > aligned.bam
  }
  output {
    File bam = "aligned.bam"
  }
  runtime {
    docker: "us.gcr.io/broad-gatk/gatk:latest"
    memory: memory
    cpu: cpu
    disks: disk
    preemptible: preemptible
  }
}

task BamToUbam {
  input {
    File bam
    File reference
    String memory
    Int cpu
    String disk
    Int preemptible
  }
  command {
    gatk RevertSam \
      -I ~{bam} \
      -O unmapped.bam \
      --SANITIZE true \
      --REMOVE_ALIGNMENT_INFORMATION true \
      --RESTORE_ORIGINAL_QUALITIES true
  }
  output {
    File ubam = "unmapped.bam"
  }
  runtime {
    docker: "us.gcr.io/broad-gatk/gatk:latest"
    memory: memory
    cpu: cpu
    disks: disk
    preemptible: preemptible
  }
}

task UbamToFastq {
  input {
    File ubam
    String prefix
    String memory
    Int cpu
    String disk
    Int preemptible
  }
  command {
    gatk SamToFastq \
      -I ~{ubam} \
      -F ~{prefix}_1.fastq \
      -F2 ~{prefix}_2.fastq \
      -FU ~{prefix}_unpaired.fastq

    gzip -c ~{prefix}_1.fastq > ~{prefix}_1.fastq.gz
    gzip -c ~{prefix}_2.fastq > ~{prefix}_2.fastq.gz
    gzip -c ~{prefix}_unpaired.fastq > ~{prefix}_unpaired.fastq.gz
  }
  output {
    File read1_fastq_gz = "~{prefix}_1.fastq.gz"
    File read2_fastq_gz = "~{prefix}_2.fastq.gz"
    File unpaired_fastq_gz = "~{prefix}_unpaired.fastq.gz"
  }
  runtime {
    docker: "us.gcr.io/broad-gatk/gatk:latest"
    memory: memory
    cpu: cpu
    disks: disk
    preemptible: preemptible
  }
  meta {
    author: "Shadi Zaheri"
  }
}