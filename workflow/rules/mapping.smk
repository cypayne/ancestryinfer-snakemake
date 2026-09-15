
# Give 4 threads to bwa and it will get 4 cores and take 
# all the memory you'd expect for those cores
# (sedna's machines are 20 core units)
# mark the bam file as temporary
rule map_reads:
    input:
        reads = [
            "results/trimmed/{sample}.1.fastq.gz",
            "results/trimmed/{sample}.2.fastq.gz" ],
        idx = rules.bwa_index.output,
    output:
        "{genome}_{sample}.bam",
    log:
        "results/logs/map_reads/{genome}_genome/{sample}.log",
    benchmark:
        "results/benchmarks/map_reads/{genome}_genome/{sample}.bmk"
    params:
        extra=get_read_group,
    threads: 4
    resources:
        mem_mb=19200,
        time="23:59:59"
    wrapper:
        "v3.13.7/bio/bwa/mem"        

