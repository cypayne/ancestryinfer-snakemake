#
# mark the sam file as temporary
rule bam_to_sam:
    input:
        "{genome}_{sample}.bam",
    output:
        temp("{genome}_{sample}.sam"),
    log:
        "results/logs/bam2sam/{genome}_{sample}.log",
    params:
        extra="-h", 
    threads: 2
    wrapper:
        "v3.13.8/bio/samtools/view"

rule rename_sam:
    input:
        sam1=genome1+"_{sample}.sam",
        sam2=genome2+"_{sample}.sam",
    output:
        renamed1="{sample}.1.fastq.gz_allchrs.par1.sam",
        renamed2="{sample}.1.fastq.gz_allchrs.par2.sam",
    shell:
        "cp {input.sam1} {output.renamed1} ; "
        "cp {input.sam2} {output.renamed2} "


# even though .parX.sam infiles aren't directly referenced 
# in the shell commands, the perl script accesses 
# these infiles
# cleans up bams and extra sams
rule samtools_to_hmm:
    input:
        unpack(get_subset_read_list),
        sam1="{sample}.1.fastq.gz_allchrs.par1.sam",
        sam2="{sample}.1.fastq.gz_allchrs.par2.sam",
        genome1=genome1,
        genome2=genome2,
    output:
        "{sample}.1.fastq.gz_allchrs.sam.hmm.combined.pass.formatted",
        "{sample}.1.fastq.gz_allchrs.sam.hmm.combinedparental.format",
    conda:
        "../envs/samtools_to_hmm.yaml",
    log:
        stderr="results/logs/variant_calling/{sample}_allchrs.stderr",
        stdout="results/logs/variant_calling/{sample}_allchrs.stdout",
    benchmark:
        "results/benchmarks/variant_calling/{sample}_allchrs.bmk",
    params:
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
        "perl {params.path}/run_samtools_to_hmm_v10.pl {input.read_list} " 
        "{input.genome1} {input.genome2} {params.read_length} "
        "{params.save_files} {params.max_align} {params.focal_chrom} "
        "{params.rec_M_per_bp} {params.path} {params.quality} > {log.stdout} "
        "2> {log.stderr} ; "

rule check_samtools_to_hmm_output:
    input: 
        expand("{sample}.1.fastq.gz_allchrs.sam.hmm.combined.pass.formatted", sample = sample_list),
        expand("{sample}.1.fastq.gz_allchrs.sam.hmm.combinedparental.format", sample = sample_list),
    output:
        "samtools_to_hmm_complete.txt", 
    shell:
        "echo 'samtools_to_hmm output file exists' >> samtools_to_hmm_complete.txt"

#        "if [ -s {input.of1} ] && [-s {input.of2} ]; then "
#        "echo '{input.of1} exists' >> samtools_to_hmm_complete.txt ;"
#        "echo '{input.of2} exists' >> samtools_to_hmm_complete.txt ;"
#        "else; echo 'at least one of the samtools_to_hmm files are empty' ;" 
#        "fi "
