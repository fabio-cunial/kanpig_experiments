version 1.0


# Resources used when run on a 16-virtual-core, 64-GB machine:
# [input, time, cost, CPU, RAM]
#
# single sample 32x CCS    
# single sample 32x ONT    
# truvari collapse 32x CCS                
# truvari collapse 32x ONT                
# bcftools merge 32x CCS   
# bcftools merge 32x ONT   
#
workflow KanpigGenotyper {
    input {
        Boolean is_singlesample
        Boolean is_male
        File input_vcf_gz
        File input_tbi
        File alignments_bam
        File alignments_bai
        File reference_fa
        File reference_fai
        Int n_cpu = 16
        Int ram_size_gb = 64
        Int sizemin = 50
        Int sizemax = 10000
        String kanpig_params_singlesample = "--chunksize 1000 --gpenalty 0.02 --hapsim 0.9999 --sizesim 0.90 --seqsim 0.85 --maxpaths 10000"
        String kanpig_params_multisample =  "--chunksize 500  --gpenalty 0.04 --hapsim 0.97"
        File ploidy_bed_female
        File ploidy_bed_male
    }
    parameter_meta {
    }
    
    call KanpigGenotyperImpl {
        input:
            is_singlesample = is_singlesample,
            is_male = is_male,
            input_vcf_gz = input_vcf_gz,
            input_tbi = input_tbi,
            alignments_bam = alignments_bam,
            alignments_bai = alignments_bai,
            reference_fa = reference_fa,
            reference_fai = reference_fai,
            n_cpu = n_cpu,
            ram_size_gb = ram_size_gb,
            sizemin = sizemin,
            sizemax = sizemax,
            kanpig_params_singlesample = kanpig_params_singlesample,
            kanpig_params_multisample = kanpig_params_multisample,
            ploidy_bed_female = ploidy_bed_female,
            ploidy_bed_male = ploidy_bed_male
    }
    
    output {
        File regenotyped_kanpig = KanpigGenotyperImpl.regenotyped_kanpig
        File regenotyped_kanpig_tbi = KanpigGenotyperImpl.regenotyped_kanpig_tbi
    }
}


task KanpigGenotyperImpl {
    input {
        Boolean is_singlesample
        Boolean is_male
        File input_vcf_gz
        File input_tbi
        File alignments_bam
        File alignments_bai
        File reference_fa
        File reference_fai
        Int n_cpu
        Int ram_size_gb
        Int sizemin
        Int sizemax
        String kanpig_params_singlesample
        String kanpig_params_multisample
        File ploidy_bed_female
        File ploidy_bed_male
    }
    parameter_meta {
    }
    
    String docker_dir = "/kanpig_experiments"
    String work_dir = "/cromwell_root/kanpig_experiments"
    String output_prefix = "kanpig_regenotyped"
    Int disk_size_gb = 200 + ceil(size(reference_fa,"GB")) + 100*ceil(size(input_vcf_gz,"GB")) + 2*ceil(size(alignments_bam,"GB"))

    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        EFFECTIVE_MEM_GB=~{ram_size_gb}
        EFFECTIVE_MEM_GB=$(( ${EFFECTIVE_MEM_GB} - 4 ))
        df -h
        
        if [ ~{is_male} == true ]; then
            PLOIDY_BED=$(echo ~{ploidy_bed_male})
        else
            PLOIDY_BED=$(echo ~{ploidy_bed_female})
        fi
        if [ ~{is_singlesample} == true ]; then
            PARAMS=$(echo ~{kanpig_params_singlesample})
        else
            PARAMS=$(echo ~{kanpig_params_multisample})
        fi
        export RUST_BACKTRACE="full"
        ${TIME_COMMAND} ~{docker_dir}/kanpig --threads $(( ${N_THREADS} - 1)) --ploidy-bed ${PLOIDY_BED} --sizemin ~{sizemin} --sizemax ~{sizemax} ${PARAMS} --reference ~{reference_fa} --input ~{input_vcf_gz} --bam ~{alignments_bam} --out tmp1.vcf.gz
        bcftools sort --max-mem ${EFFECTIVE_MEM_GB}G --output-type z tmp1.vcf.gz > ~{output_prefix}.vcf.gz
        tabix -f ~{output_prefix}.vcf.gz
    >>>

    output {
        File regenotyped_kanpig = work_dir + "/" + output_prefix + ".vcf.gz"
        File regenotyped_kanpig_tbi = work_dir + "/" + output_prefix + ".vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/kanpig_experiments"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
