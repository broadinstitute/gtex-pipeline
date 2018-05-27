task shapeit_postprocess {

    File vcf
    File vcf_index
    File phased_vcf

    String prefix
    String chrom

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

        python3 /src/shapeit_postprocess.py $chr_vcf ${phased_vcf} $PWD/${prefix}.${chrom}.phased.patched.vcf.gz
    }

    output {
        File patched_vcf = "${prefix}.${chrom}.phased.patched.vcf.gz"
        File patched_vcf_index = "${prefix}.${chrom}.phased.patched.vcf.gz.tbi"
        File log = "${prefix}.${chrom}.phased.patched.log"
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


workflow shapeit_postprocess_workflow {
    call shapeit_postprocess
}
