version 1.0

workflow GenotypeAnalysis {
    input {
	String sample
        File truth_vcf
	File truth_tbi
        File bed
	File trs
        File discovery_vcf
        File cutesv_5_vcf
        File cutesv_5_tbi
        File svjedi_5_vcf
        File svjedi_5_tbi
        File kanpig_5_vcf
        File kanpig_5_tbi
        File sniffles_5_vcf
        File sniffles_5_tbi
        File cutesv_7_vcf
        File cutesv_7_tbi
        File svjedi_7_vcf
        File svjedi_7_tbi
        File kanpig_7_vcf
        File kanpig_7_tbi
        File sniffles_7_vcf
        File sniffles_7_tbi
        File cutesv_8_vcf
        File cutesv_8_tbi
        File svjedi_8_vcf
        File svjedi_8_tbi
        File kanpig_8_vcf
        File kanpig_8_tbi
        File sniffles_8_vcf
        File sniffles_8_tbi
        File cutesv_9_vcf
        File cutesv_9_tbi
        File svjedi_9_vcf
        File svjedi_9_tbi
        File kanpig_9_vcf
        File kapig_9_tbi
        File sniffles_9_vcf
        File sniffles_9_tbi
    }

    parameter_meta {
	sample: "Name of sample. Should be found in e.g. merged VCFs"
	truth: "Dipcall SVs called from assembly"
	bed: "High confidence bed from dipcall"
	trs: "Tandem repeats bed file"
	discovery: "Sniffles discovery VCF run on sample's reads and genotyped for S5 and S7"
    }

    call GenotypeAnalysisImpl {
        input:
	    sample = sample,
	    truth = truth_vcf,
	    bed = bed,
	    trs = trs,
	    discovery = discovery_vcf,
	    cutesv_5 = cutesv_5_vcf, 
	    svjedi_5 = svjedi_5_vcf,
	    kanpig_5 = kanpig_5_vcf,
	    sniffles_5 = sniffles_5_vcf,
	    cutesv_7 = cutesv_7_vcf,
	    svjedi_7 = svjedi_7_vcf,
	    kanpig_7 = kanpig_7_vcf,
	    sniffles_7 = sniffles_7_vcf,
	    cutesv_8 = cutesv_8_vcf,
	    svjedi_8 = svjedi_8_vcf,
	    kanpig_8 = kanpig_8_vcf,
	    sniffles_8 = sniffles_8_vcf,
	    cutesv_9 = cutesv_9_vcf,
	    svjedi_9 = svjedi_9_vcf,
	    kanpig_9 = kanpig_9_vcf,
	    sniffles_9 = sniffles_9_vcf,
	    truth_tbi = truth_tbi,
	    cutesv_5_tbi = cutesv_5_tbi,
	    svjedi_5_tbi = svjedi_5_tbi,
	    kanpig_5_tbi = kanpig_5_tbi,
	    sniffles_5_tbi = sniffles_5_tbi,
	    cutesv_7_tbi = cutesv_7_tbi,
	    svjedi_7_tbi = svjedi_7_tbi,
	    kanpig_7_tbi = kanpig_7_tbi,
	    sniffles_7_tbi = sniffles_7_tbi,
	    cutesv_8_tbi = cutesv_8_tbi,
	    svjedi_8_tbi = svjedi_8_tbi,
	    kanpig_8_tbi = kanpig_8_tbi,
	    sniffles_8_tbi = sniffles_8_tbi,
	    cutesv_9_tbi = cutesv_9_tbi,
	    svjedi_9_tbi = svjedi_9_tbi,
	    kapig_9_tbi = kapig_9_tbi,
	    sniffles_9_tbi = sniffles_9_tbi,
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
	File trs
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
        File truth_tbi
        File discovery_tbi
        File cutesv_5_tbi
        File svjedi_5_tbi
        File kanpig_5_tbi
        File sniffles_5_tbi
        File cutesv_7_tbi
        File svjedi_7_tbi
        File kanpig_7_tbi
        File sniffles_7_tbi
        File cutesv_8_tbi
        File svjedi_8_tbi
        File kanpig_8_tbi
        File sniffles_8_tbi
        File cutesv_9_tbi
        File svjedi_9_tbi
        File kanpig_9_tbi
        File sniffles_9_tbi
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
	    --trs ~{trs} \
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
