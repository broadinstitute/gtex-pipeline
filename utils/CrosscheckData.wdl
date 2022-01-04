version 1.0

task CrosscheckData {
    input {
        Array[File] samples
        Array[File] samples_index
    
        File hapMap
        Int? preemptible
        Int? memoryMaybe
        String? gatkTag
    }
    String gatkTag_final = select_first([gatkTag, "4.2.4.0"])
  

    Int memoryDefault=16
    Int memoryJava=select_first([memoryMaybe,memoryDefault])
    Int memoryRam=memoryJava+2
    Int disk_size = 10 + ceil(size([hapMap], "GB"))

    File output_file="samples.crosscheck_metrics"

    parameter_meta {
        samples: {
            localization_optional: true
        }
        samples_index: {
            localization_optional: true
        }
    }

    command <<<
        set -euo pipefail

        # shellcheck disable=SC2086 # I actually want extra_arg to be split into two.
        gatk --java-options "-Xmx~{memoryJava}G" \
            CrosscheckFingerprints \
            -I ~{samples} \
            -H ~{hapMap} \
            --CALCULATE_TUMOR_AWARE_RESULTS false \
            --CROSSCHECK_BY FILE \
            --OUTPUT ~{output_file} \
    >>>
    output {
        File metrics=output_file
    }

    runtime {
            docker: "broadinstitute/gatk:" + gatkTag_final
            preemptible: select_first([preemptible, 0])
            disks: "local-disk " + disk_size + " HDD"
            bootDiskSizeGb: "16"
            memory: memoryRam + " GB"
            continueOnReturnCode: true
    }
}


workflow CrosscheckDataWF {
    input {
        Array[File] samples
        Array[File] samples_index
      
        File hapMap
        
        String? gatkTag
        
    }
    call CrosscheckData{
        input:
        gatkTag=gatkTag,
        samples=samples,
        samples_index=samples_index,
        hapMap=hapMap
    }

    output {
        File fp_metrics=CrosscheckData.metrics
    }
}