task shapeit {

    File vcf
    File vcf_index
    File pir_file
    String prefix
    String chrom

    File? sex
    File? par_bed

    Int memory
    Int disk_space
    Int num_threads
    Int num_preempt

    command {
        set -euo pipefail
        touch ${vcf_index}

        echo $(date +"[%b %d %H:%M:%S] Extracting ${chrom} from VCF")
        chr_vcf=$PWD/${chrom}.vcf.gz
        tabix -h ${vcf} ${chrom} | bgzip -c > $chr_vcf
        tabix $chr_vcf

        python3 /src/run_shapeit.py $chr_vcf ${pir_file} ${prefix} --output_dir . --num_threads ${num_threads} ${"--sex " + sex} ${"--par_bed " + par_bed}
        cat *.snp.mm > ${prefix}.${chrom}.snp.mm
        cat *.ind.mm > ${prefix}.${chrom}.ind.mm
    }

    output {
        File phased_vcf = "${prefix}.${chrom}.phased.vcf.gz"
        File snp_log = "${prefix}.${chrom}.snp.mm"
        File ind_log = "${prefix}.${chrom}.ind.mm"
    }

    runtime {
        docker: "xiaoli/shapeit2:latest"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }

    meta {
        author: "Francois Aguet"
    }
}


workflow shapeit_workflow {
    call shapeit
}

# workflow shapeit_workflow {
#
#     Array[File] chr_pir_files
#     File chr_list_file
#     Array[String] chr_list = read_lines(chr_list_file)
#
#     scatter (i in range(len(chr_list))) {
#         call shapeit { input pir_file=chr_pir_files[i], chrom=chr_list[i] }
#     }
# }
