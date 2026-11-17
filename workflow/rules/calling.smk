#
# mark the sam file as temporary
rule bam_to_sam:
    input:
        "results/mapping/{sample}/{genome}_{sample}.bam",
    output:
        temp("results/calling/{sample}/{genome}_{sample}.sam"),
    log:
        "results/logs/bam2sam/{genome}_{sample}.log",
    params:
        extra="-h", 
    threads: 2
    wrapper:
        "v3.13.8/bio/samtools/view"

rule rename_sam:
    input:
        sam1="results/calling/{sample}/"+genome1+"_{sample}.sam",
        sam2="results/calling/{sample}/"+genome2+"_{sample}.sam",
    output:
        renamed1="results/calling/{sample}/{sample}.1.fastq.gz_allchrs.par1.sam",
        renamed2="results/calling/{sample}/{sample}.1.fastq.gz_allchrs.par2.sam",
    shell:
        "cp {input.sam1} {output.renamed1} ; "
        "cp {input.sam2} {output.renamed2} "

# this step assumes you're in the same working directory as the aims file,
# so need to do some complicated funky moves to maintain organization
# even though .parX.sam infiles aren't directly referenced 
# in the shell commands, the perl script accesses these infiles
# cleans up bams and extra sams when save_files=1
rule samtools_to_hmm:
    input:
        #unpack(get_subset_read_list),
        sam1="results/calling/{sample}/{sample}.1.fastq.gz_allchrs.par1.sam",
        sam2="results/calling/{sample}/{sample}.1.fastq.gz_allchrs.par2.sam",
        genome1=f"results/ref/{genome1}",
        genome2=f"results/ref/{genome2}",
        fai=[f"results/ref/{genome1}.fai", f"results/ref/{genome2}.fai"],
        aims=f"{aims}",
        aims_mod=f"{aims}.mod",
        aims_bed=f"{aims}.mod.bed",
    output:
        "results/calling/{sample}/{sample}.1.fastq.gz_allchrs.sam.hmm.combined.pass.formatted",
        "results/calling/{sample}/{sample}.1.fastq.gz_allchrs.sam.hmm.combinedparental.format",
    conda:
        "../envs/samtools_to_hmm.yaml",
    log:
        stderr="results/logs/variant_calling/{sample}_allchrs.stderr",
        stdout="results/logs/variant_calling/{sample}_allchrs.stdout",
    benchmark:
        "results/benchmarks/variant_calling/{sample}_allchrs.bmk",
    params:
        dir="results/calling/{sample}",
        g1=lambda wc, input: os.path.abspath(input.genome1),
        g2=lambda wc, input: os.path.abspath(input.genome2),
        aims=lambda wc, input: os.path.abspath(input.aims),
        read_length=config["read_length"],
        save_files=config["retain_intermediate_files"],
        max_align=config["max_alignments"],
        focal_chrom=0,
        rec_M_per_bp=config["rec_M_per_bp"],
        path=config["program_path"],
        quality=30,
    resources:
        runtime="30h",
        mem_mb = 190000,
    threads: 1
    shell:
        "(cd {params.dir} && "
        " echo {params.aims} > current_aims_file && "
        " printf '%s\\t%s\\n' {wildcards.sample}.1.fastq.gz {wildcards.sample}.2.fastq.gz"
        "   > read_list.txt.{wildcards.sample} && "
        " perl {params.path}/run_samtools_to_hmm_v10.pl read_list.txt.{wildcards.sample} "
        "   {params.g1} {params.g2} {params.read_length} {params.save_files} "
        "   {params.max_align} {params.focal_chrom} {params.rec_M_per_bp} {params.path} {params.quality}"
        ") > {log.stdout} 2> {log.stderr}"


rule check_samtools_to_hmm_output:
    input: 
        expand("results/calling/{sample}/{sample}.1.fastq.gz_allchrs.sam.hmm.combined.pass.formatted", sample = sample_list),
        expand("results/calling/{sample}/{sample}.1.fastq.gz_allchrs.sam.hmm.combinedparental.format", sample = sample_list),
    output:
        "results/calling/samtools_to_hmm_complete.txt", 
    shell:
        "echo 'samtools_to_hmm output file exists' >> {output}"

#        "if [ -s {input.of1} ] && [-s {input.of2} ]; then "
#        "echo '{input.of1} exists' >> samtools_to_hmm_complete.txt ;"
#        "echo '{input.of2} exists' >> samtools_to_hmm_complete.txt ;"
#        "else; echo 'at least one of the samtools_to_hmm files are empty' ;" 
#        "fi "
