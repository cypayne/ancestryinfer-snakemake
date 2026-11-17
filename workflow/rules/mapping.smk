
# Give 4 threads to bwa and it will get 4 cores and take 
# all the memory you'd expect for those cores
rule map_reads:
    input:
        reads=get_reads_to_map,
        idx=rules.bwa_index.output,
    output:
        "results/mapping/{sample}/{genome}_{sample}.bam",
    log:
        "results/logs/map_reads/{genome}_genome/{sample}.log",
    #benchmark:
    #    "results/benchmarks/map_reads/{genome}_genome/{sample}.bmk"
    params:
        extra=get_read_group,
    threads: 4
    resources:
        mem_mb=19200,
        runtime="1d"
    wrapper:
        "v3.13.7/bio/bwa/mem"        

