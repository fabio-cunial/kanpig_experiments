version 1.0


#
workflow SnifflesIntersample {
    input {
        Array[File] input_snf
    }

    call JointCalling {
        input:
            input_snf = input_snf
    }

    output {
         File output_vcf_gz = JointCalling.output_vcf_gz
         File output_tbi = JointCalling.output_tbi
    }
}


task JointCalling {
    input {
        Array[File] input_snf
    }

    String docker_dir = "/kanpig_experiments"
    String work_dir = "/cromwell_root/kanpig_experiments"
    Int mem_size_gb = 32

    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}

        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        TIME_COMMAND="/usr/bin/time --verbose"
        
        INPUT_FILES=~{sep=',' input_snf}
        INPUT_FILES=$(echo ${INPUT_FILES} | tr ',' ' ')
        ${TIME_COMMAND} sniffles --input ${INPUT_FILES} --vcf tmp1.vcf
        MAX_MEM=$(( ~{mem_size_gb} - 4 ))
        bcftools sort --max-mem ${MAX_MEM}G --output-type z tmp1.vcf > joint.vcf.gz
        tabix -f joint.vcf.gz
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