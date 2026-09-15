# This turns chromosomes.tsv and scaffolds.tsv into scatter_intervals.tsv.
# Needs to have R installed and on the path.
rule dest_scatter_intervals:
    input:
        chroms=config["chromosomes"],
        scaffs=config["scaffold_groups"],
    params:
        binsize="{int_length}"
    envmodules: "R/4.0.3"
    output:
        tsv="results/scatter_config/scatters_{int_length}.tsv"
    log:
        "results/logs/dest_scatter_intervals/log_{int_length}.txt"
    script:
        "../scripts/sequence-scatter-bins.R"


