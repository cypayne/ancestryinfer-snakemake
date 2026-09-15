
## combine individuals and run ancestry_hmm
rule run_hmm:
    input:
        "samtools_to_hmm_complete.txt",
        generate_parent_list(sample_list),
        generate_hybrid_list(sample_list),
    output:
        "current.samples.list_allchrs",
        "current.samples.read.list_allchrs",
        "all.indivs.hmm.combined_allchrs.filtered_focalchroms",
    conda:
        "../envs/hmm.yaml",
    log:
        stderr="results/logs/hmm/hmm_allchrs.stderr",
        stdout="results/logs/hmm/hmm_allchrs.stdout",
    benchmark:
        "results/benchmarks/hmm/hmm_allchrs.bmk",
    params:
        minor_prior1=config["prop_genome_genome1_parent"],
        parental_counts_status=config["provide_counts"], 
        initial_admix=config["gen_initial_admix"],
        focal_chrom=0,
        read_length=config["read_length"],
        error_prior=config["per_site_error"],
        path=config["program_path"],
    resources:
        time="72:00:00",
        mem_mb = 96000,
    threads: 1
    shell:
        "perl {params.path}/combine_all_individuals_hmm_v5.pl "
        "HMM.parental.files.list_allchrs HMM.hybrid.files.list_allchrs " 
        "{params.minor_prior1} {params.parental_counts_status} {params.initial_admix} "
        "{params.focal_chrom} {params.read_length} {params.error_prior} _allchrs "
        "> {log.stdout} 2> {log.stderr} ; "


rule rchmm_to_ancestry:
    input:
        "current.samples.list_allchrs",
        "current.samples.read.list_allchrs",
    output:
        "ancestry-probs-par1_transposed_allchrs.tsv",
        "ancestry-probs-par2_transposed_allchrs.tsv",
    conda:
        "../envs/hmm.yaml",
    log:
        stderr="results/logs/hmm/rchmm2ancestry_allchrs.stderr",
        stdout="results/logs/hmm/rchmm2ancestry_allchrs.stdout",
    benchmark:
        "results/benchmarks/hmm/rchmm2ancestry_allchrs.bmk",
    params:
        save_files=0,
        focal_chrom=0,
        path=config["program_path"],
    resources:
        time="10:00:00",
        mem_mb = 190000,
    threads: 1
    shell:
        "perl {params.path}/convert_rchmm_to_ancestry_tsv_v3.pl "
        "current.samples.list current.samples.read.list " 
        "{params.save_files} {params.focal_chrom} > {log.stdout} 2> {log.stderr} "


## transpose the ancestry files
rule transpose:
    input:
        trans_ancestry1="ancestry-probs-par1_transposed_allchrs.tsv",
        trans_ancestry2="ancestry-probs-par2_transposed_allchrs.tsv",
    output:
        "ancestry-probs-par1_allchrs.tsv",
        "ancestry-probs-par2_allchrs.tsv",
    conda:
        "../envs/hmm.yaml",
    params:
        path=config["program_path"],
    shell:
        "perl {params.path}/transpose_tsv.pl {input.trans_ancestry1}; "
        "perl {params.path}/transpose_tsv.pl {input.trans_ancestry2} "
