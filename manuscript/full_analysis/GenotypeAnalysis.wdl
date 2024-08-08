version 1.0


# Resources used when run on a 16-virtual-core, 64-GB machine:
# [input, time, cost, CPU, RAM]
#
# single sample 32x CCS         3m      $ ???       660 %      1 GB 
# single sample 32x ONT         
# truvari collapse 32x CCS                
# truvari collapse 32x ONT                
# bcftools merge 32x CCS   
# bcftools merge 32x ONT   
#
workflow GenotypeAnalysis {
    input {
	String sample
        File truth
        File bed
        File discovery
        File cutesv_5
        File svjedi_5
        File kanpig_5
        File sniffles_5
        File cutesv_7
        File svjedi_7
        File kanpig_7
        File sniffles_7
        File cutesv_8
        File svjedi_8
        File kanpig_8
        File sniffles_8
        File cutesv_9
        File svjedi_9
        File kanpig_9
        File sniffles_9
    }
    parameter_meta {
    }
    
    call GenotypeAnalysisImpl {
        input:
	    sample = sample,
	    truth = truth,
	    bed = bed,
	    discovery = discovery,
	    cutesv_5 = cutesv_5, 
	    svjedi_5 = svjedi_5,
	    kanpig_5 = kanpig_5,
	    sniffles_5 = sniffles_5,
	    cutesv_7 = cutesv_7,
	    svjedi_7 = svjedi_7,
	    kanpig_7 = kanpig_7,
	    sniffles_7 = sniffles_7,
	    cutesv_8 = cutesv_8,
	    svjedi_8 = svjedi_8,
	    kanpig_8 = kanpig_8,
	    sniffles_8 = sniffles_8,
	    cutesv_9 = cutesv_9,
	    svjedi_9 = svjedi_9,
	    kanpig_9 = kanpig_9,
	    sniffles_9 = sniffles_9
    }
    
    output {
        File results_jl = GenotypeAnalysisImpl.results_jl
    }
}


task GenotypeAnalysisImpl {
    input {
	String sample
        File truth
        File bed
        File discovery
        File cutesv_5
        File svjedi_5
        File kanpig_5
        File sniffles_5
        File cutesv_7
        File svjedi_7
        File kanpig_7
        File sniffles_7
        File cutesv_8
        File svjedi_8
        File kanpig_8
        File sniffles_8
        File cutesv_9
        File svjedi_9
        File kanpig_9
        File sniffles_9
    }
    parameter_meta {
    }
    
    String docker_dir = "/kanpig_experiments"
    String work_dir = "/cromwell_root/kanpig_experiments"
    Int disk_size_gb = 24

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
        
	${TIME_COMMAND} ~{docker_dir}/analysis.py --output results.jl \
	    --sample ~{sample} \
	    --truth ~{truth} \
	    --bed ~{bed} \
	    --discovery ~{discovery} \
	    --cutesv-5 ~{cutesv_5} \
	    --svjedi-5 ~{svjedi_5} \
	    --kanpig-5 ~{kanpig_5} \
	    --sniffles-5 ~{sniffles_5} \
	    --cutesv-7 ~{cutesv_5} \
	    --svjedi-7 ~{svjedi_7} \
	    --kanpig-7 ~{kanpig_7} \
	    --sniffles-7 ~{sniffles_7} \
	    --cutesv-8 ~{cutesv_8} \
	    --svjedi-8 ~{svjedi_8} \
	    --kanpig-8 ~{kanpig_8} \
	    --sniffles-8 ~{sniffles_8} \
	    --cutesv-9 ~{cutesv_9} \
	    --svjedi-9 ~{svjedi_9} \
	    --kanpig-9 ~{kanpig_9} \
	    --sniffles-9 ~{sniffles_9}
    >>>

    output {
        File results_jl = work_dir + "/results.jl"
    }
    runtime {
        docker: "fcunial/kanpig_experiments"
        cpu: n_cpu
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
