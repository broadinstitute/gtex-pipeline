version 1.0


task ConvertPlinkToVcf{
	input {
		File bim
		File bed
		File fam

	}
	command <<<
		set -euo pipefail

		# make sure that the files all have the same basename
		# will not work if names have spaces in them
		
		bim_base=$(basename ~{bim} .bim)
		bed_base=$(basename ~{bed} .bed)
		fam_base=$(basename ~{fam} .fam)

		lines=$(echo $bim_base $bed_base $fam_base | tr ' ' '\n'| uniq | wc -l)

		if [ "${lines}" -eq "1" ]; then
			bim_bash=~{bim}
			base=${bim_bash%.bim}
			plink --bfile ${base} --recode vcf --out ~{bim_base}
		else
			echo "Error, found too many basenames in the input: $lines" 
			exit 1
		fi
	>>>

	output {
		File vcf=basename(bim,'.bim')+".vcf"
	}

	runtime {
            docker: "dnastack/plink:1.9"
            preemptible: 0
            disks: "local-disk " + ceil(size([bim,bed,fam],"GiB")+20) + " HDD"
            bootDiskSizeGb: "16"
            memory: 20 + " GB"
    }
}

workflow ConvertPlinkToVcfWF{
	call ConvertPlinkToVcf{}


    output {
        File vcf=ConvertPlinkToVcf.vcf
    }
}