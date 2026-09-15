## index the reference genomes
rule bwa_index:
    input:
        "{genome}", 
    output:
        multiext("{genome}", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        "results/logs/bwa_index_{genome}.log",
    benchmark:
        "results/benchmarks/bwa_index/bwa_index_{genome}.bmk",
    resources:
        mem_mb=36900,
    wrapper:
        "v3.13.7/bio/bwa/index"
