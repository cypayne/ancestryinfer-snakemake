## index the reference genomes
rule bwa_index:
    input:
        "results/ref/{genome}", 
    output:
        multiext("results/ref/{genome}", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        "results/logs/bwa_index/bwa_index_{genome}.log",
    benchmark:
        "results/benchmarks/bwa_index/bwa_index_{genome}.bmk",
    resources:
        mem_mb=36900,
    wrapper:
        "v3.13.7/bio/bwa/index"

## index for bcftools mpileup
rule faidx:
    input:
        "results/ref/{genome}",
    output:
        "results/ref/{genome}.fai",
#    conda:
#        "../envs/samtools.yaml",
    log:
        "results/logs/faidx/{genome}.log",
    shell:
        "samtools faidx {input} 2> {log}"
