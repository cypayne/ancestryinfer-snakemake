# ancestryinfer-snakemake

A Snakemake (≥ 8) workflow for the two-parent [ancestryinfer](https://github.com/Schumerlab/ancestryinfer) pipeline (from Schumer lab)! This pipeline infers local ancestry in hybrids from low-coverage sequencing reads, using the original ancestryinfer Perl scripts and [Ancestry_HMM](https://github.com/russcd/Ancestry_HMM).

For each sample, the workflow:

1. renames chromosomes (`_` → `-`) in both parental genomes and the AIMs file, then indexes the genomes
2. optionally trims reads with fastp 
3. maps reads to both parental genomes with `bwa mem`
4. keeps reads that map with MAPQ ≥ 30 to both genomes, then counts alleles at ancestry-informative sites (AIMs)

It then combines all samples, runs Ancestry_HMM, and writes per-site ancestry posteriors and called genotypes.

> **Run all samples together.** A site is kept for the HMM only if at least one sample in the run has reads there. The set of sites (and the results) therefore depends on which samples are in the run.

## Requirements

- Snakemake ≥ 8 and `snakemake-executor-plugin-slurm`, for example:
  `mamba create -n snakemake -c conda-forge -c bioconda snakemake snakemake-executor-plugin-slurm`
- conda/mamba. Per-rule environments are built automatically.
- A copy of the [ancestryinfer](https://github.com/Schumerlab/ancestryinfer) scripts. Set config `program_path` to that folder.
- On your `PATH` (not provided by the conda envs): `ngsutilsj`, which the ancestryinfer scripts call, and `samtools`, which the `faidx` rule calls.

## Inputs

Put these in the directory you run from:

| File | Config key | Notes |
|---|---|---|
| Parental genome FASTAs | `genome1`, `genome2` | genome1 = parent 1 (genotype `0`), genome2 = parent 2 (genotype `2`) |
| AIMs file | `provide_AIMs` | ancestryinfer format: `chrom  pos  allele1  allele2` |
| Parental counts | `provide_counts` | optional, ancestryinfer format; must have the same number of lines as the AIMs file |
| Read list | `read_list` | see below |

**Read list:** a tab-separated file with no header, one line per sample, giving the paths to read 1 and read 2:

```
reads/SAMPLE01_i7-12_read_1.fastq.gz	reads/SAMPLE01_i7-12_read_2.fastq.gz
reads/SAMPLE02_i7-40_read_1.fastq.gz	reads/SAMPLE02_i7-40_read_2.fastq.gz
```

The sample name is the text before the first `_` in the read 1 file name (`SAMPLE01` above), so every sample needs a different prefix. Only paired-end reads are supported. It's best to softlink your reads folder into the snakemake directory andreference relative paths because the pipeline uses the path to name intermediate files, and paths that are too long can throw errors.

## Configuration

Edit `config.yaml`. The main settings:

| Key | Meaning |
|---|---|
| `read_length` | length of reads, minimum distance (bp) between sites kept for the HMM |
| `prop_genome_genome1_parent` | prior proportion of the genome from parent 1 |
| `gen_initial_admix` | prior number of generations since admixture |
| `per_site_error` | per-site error rate passed to Ancestry_HMM |
| `rec_M_per_bp` | recombination rate (Morgans per bp) |
| `max_alignments` | randomly subsample each sample to this many read pairs (`0` = keep all) |
| `posterior_thresh` | minimum posterior for calling a genotype (otherwise `NA`) |
| `retain_intermediate_files` | `1` keeps BAMs, VCFs and posteriors; `0` deletes them |
| `trim_reads` | `True` trims with fastp first (options under `params: fastp:`) |

`focal_chrom_list` and `number_indiv_per_job` are not used yet. All chromosomes are run, and each sample is its own job.

## Running

Run from the repository root:

```bash
# dry run: check that inputs are found and see what will run
snakemake --configfile config.yaml -n

# on a SLURM cluster
snakemake --configfile config.yaml --profile hpcc-profiles/slurm/sherlock
```

Cluster profiles are in `hpcc-profiles/slurm/`. `sherlock` uses the SLURM executor plugin and is the one to copy for a new cluster. Set `slurm_account`, `slurm_partition`, and the default `runtime`/`mem_mb` there. Memory and runtime for individual rules are set in `workflow/rules/*.smk`.

Tips:
- **Run Snakemake itself in a small, long job** so that it isn't killed on the login node:
  ```bash
  sbatch --account=<acct> --partition=<part> --mem=4G --time=2-00:00:00 --output=driver-%j.out \
    --wrap="snakemake --configfile config.yaml --profile hpcc-profiles/slurm/sherlock"
  ```
  If it stops, submit the same command again. The run continues from the steps that finished.
- **Wrapper rules exit right away without output:** this happens when Snakemake's source cache is on a node-local disk (for example `/tmp`), because compute nodes can't see it. Point the cache at shared storage with `export XDG_CACHE_HOME=/scratch/users/$USER/.cache` before running Snakemake. If you run Snakemake in a driver job as above, put the `export` inside `--wrap`.
- **SLURM job logs** are in `.snakemake/slurm_logs/`. The rules' logs are in `results/logs/`.

## Outputs

```
results/
├── ref/         renamed genomes, bwa/fai indexes, reformatted AIMs
├── trimmed/     trimmed reads (if trim_reads: True)
├── mapping/     {sample}/ BAMs against each parent
├── calling/     {sample}/ filtered BAMs, VCFs, per-site allele counts
├── hmm/         combined HMM input, posteriors, ancestry-probs-par{1,2}_allchrs.tsv
├── genotypes/   ancestry-probs_allchrs.tsv_rec.txt
└── logs/, benchmarks/
```

- **`ancestry-probs-par{1,2}_allchrs.tsv`:** posterior probability that each site is homozygous for parent 1 or parent 2 ancestry. Samples are rows, sites (`chrom:pos`) are columns.
- **`ancestry-probs_allchrs.tsv_rec.txt`:** called genotypes. `0` = homozygous parent 1, `1` = heterozygous, `2` = homozygous parent 2, `NA` = below `posterior_thresh`.

## Citation

If you use this workflow, cite the original ancestryinfer pipeline and Ancestry_HMM. See the [ancestryinfer repository](https://github.com/Schumerlab/ancestryinfer) for references.
