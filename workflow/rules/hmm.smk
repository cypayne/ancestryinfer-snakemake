
## combine individuals and run ancestry_hmm
rule run_hmm:
    input:
        hyb=expand("results/calling/{sample}/{sample}.1.fastq.gz_allchrs.sam.hmm.combined.pass.formatted", sample=sample_list),
        par=expand("results/calling/{sample}/{sample}.1.fastq.gz_allchrs.sam.hmm.combinedparental.format", sample=sample_list),
        aims=f"{aims}",
        counts=[config["provide_counts"]] if config["provide_counts"] else [],
    output:
        "results/hmm/current.samples.list_allchrs",
        "results/hmm/current.samples.read.list_allchrs",
        "results/hmm/all.indivs.hmm.combined_allchrs.filtered_focalchroms",
    conda:
        "../envs/hmm.yaml",
    log:
        stderr="results/logs/hmm/hmm_allchrs.stderr",
        stdout="results/logs/hmm/hmm_allchrs.stdout",
    benchmark:
        "results/benchmarks/hmm/hmm_allchrs.bmk",
    params:
        hyb=lambda wc, input: [os.path.abspath(f) for f in input.hyb],
        par=lambda wc, input: [os.path.abspath(f) for f in input.par],
        aims=lambda wc, input: os.path.abspath(input.aims),
        counts=os.path.abspath(config["provide_counts"]) if config["provide_counts"] else 0,
        minor_prior1=config["prop_genome_genome1_parent"],
        parental_counts_status=config["provide_counts"], 
        initial_admix=config["gen_initial_admix"],
        focal_chrom=0,
        read_length=config["read_length"],
        error_prior=config["per_site_error"],
        path=config["program_path"],
    resources:
        runtime="72h",
        mem_mb = 250000,
    threads: 1
    shell:
        "(cd results/hmm && "
        " echo {params.aims} > current_aims_file && "
        " printf '%s\\n' {params.hyb} > HMM.hybrid.files.list_allchrs && "
        " printf '%s\\n' {params.par} > HMM.parental.files.list_allchrs && "
        " perl {params.path}/combine_all_individuals_hmm_v5.pl "
        "   HMM.parental.files.list_allchrs HMM.hybrid.files.list_allchrs "
        "   {params.minor_prior1} {params.counts} {params.initial_admix} 0 "
        "   {params.read_length} {params.error_prior} _allchrs"
        ") > {log.stdout} 2> {log.stderr}"


rule rchmm_to_ancestry:
    input:
        "results/hmm/current.samples.list_allchrs",
        "results/hmm/current.samples.read.list_allchrs",
    output:
        "results/hmm/ancestry-probs-par1_transposed_allchrs.tsv",
        "results/hmm/ancestry-probs-par2_transposed_allchrs.tsv",
    conda:
        "../envs/hmm.yaml",
    log:
        stderr="results/logs/hmm/rchmm2ancestry_allchrs.stderr",
        stdout="results/logs/hmm/rchmm2ancestry_allchrs.stdout",
    benchmark:
        "results/benchmarks/hmm/rchmm2ancestry_allchrs.bmk",
    params:
        save_files=config["retain_intermediate_files"],
        focal_chrom=0,
        path=config["program_path"],
    resources:
        runtime="10h",
        mem_mb = 190000,
    threads: 1
    shell:
        "(cd results/hmm && "
        " perl {params.path}/convert_rchmm_to_ancestry_tsv_v3.pl "
        "   current.samples.list current.samples.read.list "
        "   {params.save_files} {params.focal_chrom}"
        ") > {log.stdout} 2> {log.stderr}"

## transpose the ancestry files
rule transpose:
    input:
        trans_ancestry1="results/hmm/ancestry-probs-par1_transposed_allchrs.tsv",
        trans_ancestry2="results/hmm/ancestry-probs-par2_transposed_allchrs.tsv",
    output:
        "results/hmm/ancestry-probs-par1_allchrs.tsv",
        "results/hmm/ancestry-probs-par2_allchrs.tsv",
    conda:
        "../envs/hmm.yaml",
    params:
        path=config["program_path"],
    shell:
        "perl {params.path}/transpose_tsv.pl {input.trans_ancestry1}; "
        "perl {params.path}/transpose_tsv.pl {input.trans_ancestry2} "


## calculate ancestry proportions (hybrid indices) for each sample
rule calc_hybrid_indices:
    input:
        ancestry1="results/hmm/ancestry-probs-par1_allchrs.tsv",
        ancestry2="results/hmm/ancestry-probs-par2_allchrs.tsv",
    output:
        hi_file="results/hmm/hybrid_index_file",
    params:
        path=config["program_path"],
    shell:
        "perl {params.path}/parsetsv_ancestry_v2.pl {input.ancestry1} {input.ancestry2} > {output.hi_file} " 



## genotype ancestry probabilities 
## homozygous par1 coded as 0, homozygous par2 coded as 2
rule genotype:
    input:
        ancestry1="results/hmm/ancestry-probs-par1_allchrs.tsv",
        ancestry2="results/hmm/ancestry-probs-par2_allchrs.tsv",
    output:
        geno_file="results/genotypes/ancestry-probs_allchrs.tsv_rec.txt",
    conda:
        "../envs/hmm.yaml",
    log:
        stderr="results/logs/genotypes/genotype.stderr",
        stdout="results/logs/genotypes/genotype.stdout",
    params:
        path=config["program_path"],
        posterior_thresh=config["posterior_thresh"],
    shell:
        "perl {params.path}/parsetsv_to_genotypes_v2.pl {input.ancestry1} "
        "{input.ancestry2} {params.posterior_thresh} {output.geno_file} > "
        "{log.stdout} 2> {log.stderr} "


## quantify missingness and filter out samples with too much missing data
rule missing_calls:
    input:
        geno_file="results/genotypes/ancestry-probs_allchrs.tsv_rec.txt",
    output:
        miss_counts_file="results/genotypes/missing_call_counts.txt",
        filtered_geno_file="results/genotypes/pass-miss-threshold.ancestry-probs_allchrs.tsv_rec.txt",
    params:
        max_pct=config["min_percent_missing_calls"],
    shell:
        r"""
        awk -F'\t' \
            -v thr={params.max_pct} \
            -v counts_out={output.miss_counts_file} \
            -v filt_out={output.filtered_geno_file} '
            BEGIN {{print "sample\tn_missing\tpct_missing" > counts_out}}
            NR == 1 {{print > filt_out; next}}
            {{
                missing = 0
                total = NF - 1
                for (i = 2; i <= NF; i++) {{ if ($i == "NA") missing++ }}
                pct = (total > 0) ? 100 * missing / total : 0
                printf "%s\t%d\t%.2f\n", $1, missing, pct > counts_out
                if (pct < thr) print > filt_out
            }}' \
            {input.geno_file}
        """


## cleanup everything in results/calling/{sample} that isn't necessary to rerun hmm steps  
rule cleanup:
    input:
        expand("results/calling/{sample}/{sample}.1.fastq.gz_allchrs.sam.hmm.combined.pass.formatted", sample=sample_list),
        expand("results/calling/{sample}/{sample}.1.fastq.gz_allchrs.sam.hmm.combinedparental.format", sample=sample_list),
    output:
        "results/calling/cleanup_complete.txt",
    log:
        "results/logs/cleanup/cleanup.log",
    params:
        dirs=expand("results/calling/{sample}", sample=sample_list),
    shell:
        r"""
        for d in {params.dirs}; do
            find "$d" -mindepth 1 -type f \
                ! -name '*.sam.hmm.combined.pass.formatted' \
                ! -name '*.sam.hmm.combinedparental.format' \
                -print -delete
            find "$d" -mindepth 1 -type d -empty -print -delete
        done > {log} 2>&1
        echo "cleanup finished $(date)" > {output}
        """
