
#
rule trim_reads_pe:
    input:
        unpack(get_fastq),
    output:
        read1="results/trimmed/{sample}.1.fastq.gz",
        read2="results/trimmed/{sample}.2.fastq.gz",
        html="results/qc/fastp/{sample}.html",
        json="results/qc/fastp/{sample}.json",
    conda:
        "../envs/fastp.yaml",
    log:
        out="results/logs/trimmed_pe/{sample}.log",
        err="results/logs/trimmed_pe/{sample}.err",
    benchmark:
        "results/benchmarks/trimmed_pe/{sample}.bmk",
    params:
        trim_settings=config["params"]["fastp"],
    shell:
        " fastp -i {input.read1} -I {input.read2} "
        "       -o {output.read1} -O {output.read2} "
        "       -h {output.html} -j {output.json} "
        "  {params.trim_settings} > {log.out} 2> {log.err} "
