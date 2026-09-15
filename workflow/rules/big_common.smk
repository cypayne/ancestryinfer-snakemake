import os
import warnings
import pandas as pd
from snakemake.utils import validate
from snakemake.utils import min_version

# The two lines below can be used to get the config information
# for testing.  Change config file path as appropriate
#  from snakemake.io import _load_configfile
# config = _load_configfile("config.yaml")

min_version("5.18.0")
container: "continuumio/miniconda3:4.8.2"

###### Config file and sample sheets #####

read_list = pd.read_table(config["read_list"], dtype=str).set_index(
    ["read1", "read2"], drop=False
)
read_list.index = read_list.index.set_levels(
    [i.astype(str) for i in read_list.index.levels]
)  # enforce str in index

### Eric's addition to genotype over the chromosomes and scaffold_groups
chromosomes = pd.read_table(config["chromosomes"]).set_index("chrom", drop=False)

scaffold_groups = pd.read_table(config["scaffold_groups"]).set_index("id", drop=False)

# column 1 of the scaffold group file must be "id" and
# column 2 must be "chrom". These headers are used to downstream. 
scaff_cols = list(scaffold_groups.columns)
if scaff_cols[0] != 'id' or scaff_cols[1] != 'chrom': 
    raise Exception("Column order is important in the scaffold_groups file. The first column must be 'id' and the second column must be 'chrom'.")

# get a list of just the unique values of the scaffold_group and of the chromosomes
unique_scaff_groups = list(scaffold_groups.id.unique())
unique_chromosomes = list(chromosomes.chrom.unique())  # don't need to unique it, but I do anyway

# finally, get all the scatter groups, indexed two different ways.
scatter_wc_constraint="scat_0[0-9]*"  # this is just here for the case where `scatter_intervals_file: ""`
if config["scatter_intervals_file"] != "":
    scatter_groups = pd.read_table(config["scatter_intervals_file"]).set_index("id", drop=False)
    validate(scatter_groups, schema="../schemas/scatter_intervals.schema.yaml")
    scatter_cols = list(scatter_groups.columns)
    if scatter_cols[0] != 'id' or scatter_cols[1] != 'scatter_idx' or scatter_cols[2] != 'chrom' or scatter_cols[3] != 'start' or scatter_cols[4] != 'end' or scatter_cols[5] != 'scatter_length':
        raise Exception("Column order is important in the scaffold_groups file.  The columns must be in order: id, scatter_idx, chrom, start, end, scatter_length.")
    unique_scats = list(scatter_groups.scatter_idx.unique())
    scatter_wc_constraint="|".join(unique_scats)


##### Wildcard constraints #####
# we have to deal with the cases where scaff groups or chroms might be empty
if not unique_scaff_groups:
    sg_constraint = "CrazyScaffGroupName"  # this is sort of a hack.  Just put something that won't be a chromosome
else:
    sg_constraint="|".join(unique_scaff_groups)

if not unique_chromosomes:
    chrom_constraint="CrazyChromName"
else:
    chrom_constraint="|".join(unique_chromosomes)

wildcard_constraints:
    genome="|".join(genome_list),
    sample="|".join(sample_list),
    read_list="|".join(read_list),
    chromo=chrom_constraint,
    scaff_group=sg_constraint,
    sg_or_chrom="|".join(unique_scaff_groups + unique_chromosomes),
    scatter=scatter_wc_constraint


##### Helper functions #####

def get_genome_paths(wildcards):
    """grab genome files""" 
    genomes = {genome_list}
    if len(genomes) == 2:
        return {"genome1": genomes[0], "genome2": genomes[1]}
    return {"genome1": genomes[0]}

def get_fastq(wildcards):
    """Get fastq files of given sample-unit."""
    fastqs = read_list.loc[["read1", "read2"]].dropna()
    if len(fastqs) == 2:
        return {"read1": fastqs.read1, "read2": fastqs.read2}
    return {"read1": fastqs.read1}

def is_single_end(wildcards):
    """Return True if sample-unit is single end."""
    if config["read_type"] == SE: return False
    else return True 

def get_read_group(wildcards):
    """Denote sample name and platform in read group."""
    return r"-R '@RG\tID:{sample}_{sample_id}_{library}_{flowcell}_{lane}_{barcode}\tSM:{sample_id}\tPL:{platform}\tLB:{library}\tPU:{flowcell}.{lane}.{barcode}'".format(
        sample=wildcards.sample,
        sample_id=units.loc[(wildcards.sample, wildcards.unit), "sample_id"],
        platform=units.loc[(wildcards.sample, wildcards.unit), "platform"],
        library=units.loc[(wildcards.sample, wildcards.unit), "library"],
        flowcell=units.loc[(wildcards.sample, wildcards.unit), "flowcell"],
        lane=units.loc[(wildcards.sample, wildcards.unit), "lane"],
        barcode=units.loc[(wildcards.sample, wildcards.unit), "barcode"],
    )

def get_trimmed_reads(wildcards):
    """Get trimmed reads of given sample-unit."""
    if not is_single_end(**wildcards):
        # paired-end sample
        return expand(
            "results/bqsr-round-{{bqsr_round}}/trimmed/{sample}---{unit}.{group}.fastq.gz", group=[1, 2], **wildcards
        )
    # single end sample
    return "results/bqsr-round-{{bqsr_round}}/trimmed/{sample}---{unit}.fastq.gz".format(**wildcards)

def get_sample_bams(wildcards):
    """Get all aligned reads of given sample."""
    return expand(
        "results/bqsr-round-{{bqsr_round}}/mkdup/{sample}---{unit}.bam",
        sample=wildcards.sample,
        unit=units.loc[wildcards.sample].unit,
    )

def get_regions_param(regions=config["processing"].get("restrict-regions"), default=""):
    if regions:
        params = "--intervals '{}' ".format(regions)
        padding = config["processing"].get("region-padding")
        if padding:
            params += "--interval-padding {}".format(padding)
        return params
    return default


## ancestryHMM additions

genome_list = list(config["genome1"],config["genome2"]) 
