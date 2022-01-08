version 1.0

task CrosscheckData {
    input {
        Array[File] samples
        Array[File] samples_index
    
        File hapMap
        Int? preemptible
        Int? memoryMaybe
        String? gatkTag
        Int threads
    }
    String gatkTag_final = select_first([gatkTag, "4.2.4.0"])
  

    Int memoryDefault=16
    Int memoryJava=select_first([memoryMaybe,memoryDefault])
    Int memoryRam=memoryJava+2
    Int disk_size = 10 + ceil(size([hapMap], "GB"))

    String output_name="samples.crosscheck_metrics"

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

        gatk --java-options "-Xmx~{memoryJava}G" \
            CrosscheckFingerprints \
            -I ~{sep=" -I " samples} \
            -NUM_THREADS ~{threads} \
            -H ~{hapMap} \
            --CALCULATE_TUMOR_AWARE_RESULTS false \
            --CROSSCHECK_BY FILE \
            --OUTPUT ~{output_name} \
    >>>
    output {
        File metrics=output_name
    }

    runtime {
            docker: "broadinstitute/gatk:" + gatkTag_final
            preemptible: select_first([preemptible, 0])
            disks: "local-disk " + disk_size + " HDD"
            bootDiskSizeGb: "16"
            memory: memoryRam + " GB"
            cpu: 2
            continueOnReturnCode: true
    }
}


workflow CrosscheckDataWF {
    input {
        Array[File] samples
        Array[File] samples_index
      
        File hapMap
        String? gatkTag
        Int? threads
    }
    Int length_min_forty = if length(samples) >= 40 then 40 else length(samples) 

    Int threads_final = select_first([threads,length_min_forty])

    call CrosscheckData{
        input:
        gatkTag=gatkTag,
        samples=samples,
        threads=threads_final,
        samples_index=samples_index,
        hapMap=hapMap
    }

    output {
        File fp_metrics=CrosscheckData.metrics
    }
}