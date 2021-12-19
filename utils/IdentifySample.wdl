version 1.0

import "PicardMetrics.wdl" as p

task IdentifySample {
    input {
        File sample
        File sample_index
        File vcf
        File vcf_index

        String? expected_sample_name

        File hapMap
        Int? preemptible
        Int? memoryMaybe
        String? gatkTag
    }
    String gatkTag_final = select_first([gatkTag, "4.2.4.0"])
    String expected_sample_name_final = select_first([expected_sample_name, ""])


    Int memoryDefault=16
    Int memoryJava=select_first([memoryMaybe,memoryDefault])
    Int memoryRam=memoryJava+2
    Int disk_size = 10 + ceil(size([hapMap, vcf], "GB"))

    parameter_meta {
        sample: {
            localization_optional: true
        }
        sample_index: {
            localization_optional: true
        }
        vcf: {
            localization_optional: true
        }
        vcf_index: {
            localization_optional: true
        }
    }

    command <<<
        set -euo pipefail

        if [ "~{expected_sample_name_final}" != "" ]; then
            printf '~{sample}\t~{expected_sample_name_final}\n' > sample_map.txt
            extra_arg="--INPUT_SAMPLE_MAP sample_map.txt"
        else
            extra_arg=""    
        fi

        gatk --java-options "-Xmx~{memoryJava}G" \
            CrosscheckFingerprints \
            -I ~{sample} \
            -SI ~{vcf} \
            "${extra_arg}" \
            -H ~{hapMap} \
            --CALCULATE_TUMOR_AWARE_RESULTS false \
            --CROSSCHECK_MODE CHECK_ALL_OTHERS \
            --CROSSCHECK_BY SAMPLE \
            --OUTPUT sample.crosscheck_metrics \
    >>>
    output {
        File metrics="sample.crosscheck_metrics"
    }

    runtime {
            docker: "broadinstitute/gatk:" + gatkTag_final
            preemptible: select_first([preemptible, 0])
            disks: "local-disk " + disk_size + " HDD"
            bootDiskSizeGb: "16"
            memory: memoryRam + " GB"
    }
}

task ClusterMetrics {
    input {
        File fp_metrics
        String? gatkTag
    }
    String gatkTag_final = select_first([gatkTag, "4.2.4.0"])


    Int memoryJava=16
    Int memoryRam=memoryJava+2
    Int disk_size = 15 

    
    command <<<
        set -euo pipefail

        gatk --java-options "-Xmx~{memoryJava}G" \
            ClusterFingerprintMetrics \
            -I ~{fp_metrics} \
            --OUTPUT sample.clustered.crosscheck_metrics \
    >>>
    output {
        File metrics="sample.crosscheck_metrics"
    }

    runtime {
            docker: "broadinstitute/gatk:" + gatkTag_final
            disks: "local-disk " + disk_size + " HDD"
            bootDiskSizeGb: "16"
            memory: memoryRam + " GB"
            continueOnReturnCode: [0,3]
    }
}

workflow IdentifySampleWF{
    input {
        String? gatkTag
    }
    call IdentifySample{
        input:
        gatkTag=gatkTag
    }

    call ClusterMetrics{input:
        gatkTag=gatkTag,
        fp_metrics=IdentifySample.metrics
    }

    call p.LoadPicardMetricsWF as picard_metric{
        input:
        picard_metrics=ClusterMetrics.metrics,
        requested_metric="RIGHT_GROUP_VALUE",
        requested_row=0
    }
    output {
        File fp_metrics=IdentifySample.metrics
        File fp_clustered=ClusterMetrics.metrics
        String match_group=picard_metric.requested_value
    }
}