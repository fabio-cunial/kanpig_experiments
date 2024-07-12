version 1.0


# Resources used when run on a 8-virtual-core, 32-GB machine:
# [input, time, cost, CPU, RAM]
#
# 32x CCS   h   $      %     GB
# 32x ONT   h   $      %     GB
#
workflow LrcallerGenotyper {
    input {
        File input_vcf_gz
        File input_tbi
        File alignments_bam
        File alignments_bai
        File reference_fa
        File reference_fai
        String genotyper = "joint"
    }
    parameter_meta {
        genotyper: "Options: joint, ad, va, multi. Default: joint."
    }
    
    call LrcallerGenotyperImpl {
        input:
            input_vcf_gz = input_vcf_gz,
            input_tbi = input_tbi,
            alignments_bam = alignments_bam,
            alignments_bai = alignments_bai,
            reference_fa = reference_fa,
            reference_fai = reference_fai,
            genotyper = genotyper
    }
    
    output {
        #File regenotyped_lrcaller = LrcallerGenotyperImpl.regenotyped_lrcaller
        #File regenotyped_lrcaller_tbi = LrcallerGenotyperImpl.regenotyped_lrcaller_tbi
        File monitor_log = LrcallerGenotyperImpl.monitor_log
    }
}


task LrcallerGenotyperImpl {
    input {
        File input_vcf_gz
        File input_tbi
        File alignments_bam
        File alignments_bai
        File reference_fa
        File reference_fai
        String genotyper
    }
    parameter_meta {
    }
    
    String docker_dir = "/kanpig_experiments"
    String work_dir = "/cromwell_root/kanpig_experiments"
    String output_prefix = "lrcaller_regenotyped"
    Int disk_size_gb = 500 + ceil(size(reference_fa,"GB")) + 10*ceil(size(input_vcf_gz,"GB")) + ceil(size(alignments_bam,"GB"))

    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        df -h
        
        # Starting resource monitoring
        export MONITOR_MOUNT_POINT=~{work_dir}
        MONITOR_FILE="./monitor.log"
        bash ~{docker_dir}/vm_local_monitoring_script.sh &> ${MONITOR_FILE} &
        MONITOR_JOB=$(ps -aux | grep -F 'vm_local_monitoring_script.sh' | head -1 | awk '{print $2}')

        ${TIME_COMMAND} lrcaller --number_of_threads ${N_THREADS} --dyn-w-size --genotyper ~{genotyper} --fa ~{reference_fa} ~{alignments_bam} ~{input_vcf_gz} ~{output_prefix}.vcf 2> /dev/null && echo 0 || echo 1
        
        # Stopping resource monitoring
        kill ${MONITOR_JOB}
        tail -n 100 ${MONITOR_FILE}
        df -h
        ls -laht
        
        # Outputting
        bgzip ~{output_prefix}.vcf
        tabix -f ~{output_prefix}.vcf.gz
    >>>

    output {
        #File regenotyped_lrcaller = work_dir + "/" + output_prefix + ".vcf.gz"
        #File regenotyped_lrcaller_tbi = work_dir + "/" + output_prefix + ".vcf.gz.tbi"
        File monitor_log = work_dir + "/monitor.log"
    }
    runtime {
        docker: "fcunial/kanpig_experiments"
        cpu: 8
        memory: "128GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
