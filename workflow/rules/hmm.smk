
## combine individuals and run ancestry_hmm
rule run_hmm:
    input:
        hyb=expand("results/calling/{sample}/{sample}.1.fastq.gz_allchrs.sam.hmm.combined.pass.formatted", sample=sample_list),
        par=expand("results/calling/{sample}/{sample}.1.fastq.gz_allchrs.sam.hmm.combinedparental.format", sample=sample_list),
        aims=f"{aims}",
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
        mem_mb = 96000,
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
