version 1.0


#
workflow SnifflesSingleSample {
    input {
        File input_bam
        File input_bai
        String input_id
        File tandems_bed
        File reference_fa
    }

    call SingleSampleCalling {
        input:
            input_bam = input_bam,
            input_bai = input_bai,
            sample_id = input_ids,
            tandems_bed = tandems_bed,
            reference_fa = reference_fa
    }

    output {
         File output_vcf_gz = SingleSampleCalling.output_vcf_gz
         File output_tbi = SingleSampleCalling.output_tbi
    }
}


task SingleSampleCalling {
    input {
        File input_bam
        File input_bai
        String sample_id
        File tandems_bed
        File reference_fa
    }
    
    String docker_dir = "/hgsvc2"
    String work_dir = "/cromwell_root/hgsvc2"

    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
    
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        TIME_COMMAND="/usr/bin/time --verbose"

        ${TIME_COMMAND} sniffles --threads ${N_THREADS} \
                 --input ~{input_bam} \
                 --tandem-repeats ~{tandems_bed} \
                 --reference ~{reference_fa} \
                 --sample-id ~{sample_id} \
                 --vcf ~{sample_id}.sniffles.vcf \
                 --snf ~{sample_id}.sniffles.snf
        ls -laht
        tree
    >>>

    output {
        File snf = work_dir + "/" + sample_id + ".sniffles.snf"
        File vcf = work_dir + "/" + sample_id + ".sniffles.vcf"
    }
    runtime {
        docker: "fcunial/callset_integration"
        cpu: 16
        memory: "32GB"  # Arbitrary
        disks: "local-disk 50 HDD"  # Arbitrary
        preemptible: 0
    }
}


task JointCalling {
    input {
        Array[File] input_snf
    }

    String docker_dir = "/hgsvc2"
    String work_dir = "/cromwell_root/hgsvc2"

    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}

        GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
        GSUTIL_DELAY_S="600"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        TIME_COMMAND="/usr/bin/time --verbose"
        
        INPUT_FILES=~{sep=',' input_snf}
        INPUT_FILES=$(echo ${INPUT_FILES} | tr ',' ' ')
        ${TIME_COMMAND} sniffles --input ${INPUT_FILES} --vcf joint.vcf
        bgzip joint.vcf
        tabix joint.vcf.gz
        ls -laht
        tree
    >>>

    output {
        File output_vcf_gz = work_dir + "/joint.vcf.gz"
        File output_tbi = work_dir + "/joint.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/kanpig_experiments"
        cpu: 16
        memory: "32GB"  # Arbitrary
        disks: "local-disk 50 HDD"  # Arbitrary
        preemptible: 0
    }
}