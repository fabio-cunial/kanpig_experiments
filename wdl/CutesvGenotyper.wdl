version 1.0


# Resources used when run on a 16-virtual-core, 64-GB machine:
# [input, time, cost, CPU, RAM]
#
# single sample 32x CCS     1h      $ 1.48     33 %       2 GB
# single sample 32x ONT     3h      $ 2.04    160 %       3 GB
# truvari collapse 32x CCS                
# truvari collapse 32x ONT                
# bcftools merge 32x CCS    1h      $ 1       170 %      10 GB
# bcftools merge 32x ONT    1h      $ 1.50    500 %      10 GB
#
workflow CutesvGenotyper {
    input {
        File input_vcf_gz
        File input_tbi
        File alignments_bam
        File alignments_bai
        File reference_fa
        File reference_fai
        String cutesv_params_ccs = "--max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5 --merge_ins_threshold 500 --merge_del_threshold 500"
        String cutesv_params_ont = "--max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 --merge_ins_threshold 500 --merge_del_threshold 500"
        Int is_ont
        Int cutesv_minsupport = 10
        Int svlen_min = 30
        Int svlen_max = -1
    }
    parameter_meta {
        cutesv_params_ccs: "Default values are taken from the command-line help."
        cutesv_params_ont: "Default values are taken from the command-line help."
        cutesv_minsupport: "Default=10. From the source code it seems to be used only in discovery, so it is probably irrelevant for genotyping."
        svlen_min: "Default=30"
        is_ont: "0/1"
    }
    
    call CutesvGenotyperImpl {
        input:
            input_vcf_gz = input_vcf_gz,
            input_tbi = input_tbi,
            alignments_bam = alignments_bam,
            alignments_bai = alignments_bai,
            reference_fa = reference_fa,
            reference_fai = reference_fai,
            cutesv_params_ccs = cutesv_params_ccs,
            cutesv_params_ont = cutesv_params_ont,
            is_ont = is_ont,
            cutesv_minsupport = cutesv_minsupport,
            svlen_min = svlen_min,
            svlen_max = svlen_max
    }
    
    output {
        File regenotyped_cutesv = CutesvGenotyperImpl.regenotyped_cutesv
        File regenotyped_cutesv_tbi = CutesvGenotyperImpl.regenotyped_cutesv_tbi
    }
}


task CutesvGenotyperImpl {
    input {
        File input_vcf_gz
        File input_tbi
        File alignments_bam
        File alignments_bai
        File reference_fa
        File reference_fai
        String cutesv_params_ccs
        String cutesv_params_ont
        Int is_ont
        Int cutesv_minsupport
        Int svlen_min
        Int svlen_max
    }
    parameter_meta {
    }
    
    String docker_dir = "/kanpig_experiments"
    String work_dir = "/cromwell_root/kanpig_experiments"
    String output_prefix = "cutesv_regenotyped"
    Int disk_size_gb = 200 + ceil(size(reference_fa,"GB")) + 100*ceil(size(input_vcf_gz,"GB")) + 2*ceil(size(alignments_bam,"GB"))

    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        df -h
        
        if [ ~{is_ont} -eq 1 ]; then
            PARAMS=$(echo ~{cutesv_params_ont})
        else
            PARAMS=$(echo ~{cutesv_params_ccs})
        fi
        mkdir ./cutesv_tmp
        ${TIME_COMMAND} cuteSV --threads ${N_THREADS} -Ivcf ~{input_vcf_gz} ${PARAMS} --genotype --min_support ~{cutesv_minsupport} --min_size ~{svlen_min} --max_size ~{svlen_max} ~{alignments_bam} ~{reference_fa} ~{output_prefix}.vcf ./cutesv_tmp
        rm -rf ./cutesv_tmp
        bgzip ~{output_prefix}.vcf
        tabix -f ~{output_prefix}.vcf.gz
    >>>

    output {
        File regenotyped_cutesv = work_dir + "/" + output_prefix + ".vcf.gz"
        File regenotyped_cutesv_tbi = work_dir + "/" + output_prefix + ".vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/kanpig_experiments"
        cpu: 2
        memory: "64GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
