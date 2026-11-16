# 

import os
import warnings
import pandas as pd
from snakemake.utils import validate
from snakemake.utils import min_version

min_version("8.0")
container: "continuumio/miniconda3:4.8.2"

###### Config file and sample sheets #####

# format the read list and add a column called sample, which is the first
# value when the read name is split by "_" 
read_list = pd.read_table(config["read_list"], dtype=str, names=["read1", "read2"])
read_list['sample'] = read_list['read1'].str.split('/').str[-1].str.split('_').str[0]

read_list = read_list.set_index(["sample"],drop=False)

sample_list = list(read_list["sample"].unique())
print(sample_list, file=sys.stderr)


##### Wildcard constraints #####

wildcard_constraints:
    sample="|".join(sample_list),
    read_list="|".join(read_list)


##### Helper functions #####

def get_fastq(wildcards):
    """Get fastq files of given sample-unit."""
    fastqs = read_list.loc[(wildcards.sample), ["read1", "read2"]].dropna()
    return {"read1": fastqs.read1, "read2": fastqs.read2}

def get_reads_to_map(wildcards):
    """Trimmed reads if config trim_reads is True, otherwise raw reads."""
    if config.get("trim_reads", True):
        return expand("results/trimmed/{sample}.{r}.fastq.gz",
                      sample=wildcards.sample, r=[1, 2])
    fq = get_fastq(wildcards)
    return [fq["read1"], fq["read2"]]

def get_read_group(wildcards):
    """Denote sample name and platform in read group."""
    return r"-M -R '@RG\tID:{sample}_{sample}_lib_flowcell_lane_barcode\tSM:{sample}\tPL:illumina\tLB:hyblib1\tPU:flowcell.lane.barcode'".format(
        sample=wildcards.sample,
    )

def get_subset_read_list(wildcards):
    """Return a list of id files, each containing one sample: """
    """trimmed_fq1_path\ttrimmed_fq2_path"""
    with open("read_list.txt."+wildcards.sample, 'w') as f:
        f.write(wildcards.sample+".1.fastq.gz\t"+wildcards.sample+".2.fastq.gz")
    return {"read_list": expand("read_list.txt.{sample}", sample=wildcards.sample)}

def get_sample_bams(wildcards):
    """Get all aligned reads of given sample."""
    return expand(
        "results/bqsr-round-{{bqsr_round}}/mkdup/{sample}---{unit}.bam",
        sample=wildcards.sample,
        unit=units.loc[wildcards.sample].unit,
    )

## write list of parental files to a hybrid list file
def generate_parent_list(samples):
    with open("HMM.parental.files.list_allchrs", 'w') as f:
        for s in samples:
            f.write(s+".1.fastq.gz_allchrs.sam.hmm.combinedparental.format\n")
        return "HMM.parental.files.list_allchrs" 

## write list of hybrid files to a hybrid list file
def generate_hybrid_list(samples):
    with open("HMM.hybrid.files.list_allchrs", 'w') as f:
        for s in samples:
            f.write(s+".1.fastq.gz_allchrs.sam.hmm.combined.pass.formatted\n")
    return "HMM.hybrid.files.list_allchrs" 


