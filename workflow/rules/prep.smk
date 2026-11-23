
# check read count per sample before anything
# this performs a rough check that fastqs have
# a minimum number of reads, outputs new read_list
# with files that meet this requirement
rule min_read_count_check:
    input:
        read_list = config["read_list"],
    output:
        read_counts = "results/read_counts/read_counts.1.txt",
        new_read_list = "results/read_counts/new_read_list.txt",
    params:
        min_rc = config["min_read_count"],
    shell:
        r"""
        awk '{{print $1}}' {input.read_list} | while IFS= read -r file; do
            ls -lL "$file" | awk '{{printf("%s\t%s\n", $NF, $5/1000.0)}}'
        done > {output.read_counts}
        awk -F'\t' -v min={params.min_rc} 'NR==FNR {{ if ($2 > min) keep[$1] = 1; next }} ($1 in keep)' \
            {output.read_counts} {input.read_list} > {output.new_read_list}
        """


# rename chromosome names to replace underscores with dashes
# unfortunately cannot have underscores in chrom names
rule rename_chrom_names:
    input: 
        og_g1 = og_genome1,
        og_g2 = og_genome2,
    output: 
        renamed_g1 = f"results/ref/{og_genome1}_chr-renamed",
        renamed_g2 = f"results/ref/{og_genome2}_chr-renamed",
    shell:
        "perl -p -e 's/_/-/g' {input.og_g1} > {output.renamed_g1} ; "
        "perl -p -e 's/_/-/g' {input.og_g2} > {output.renamed_g2} "


# reformat the AIMs file for downstream compatibility
# make "current_aims_file"
rule reformat_aims:
    input:
        provided_aims=config["provide_AIMs"],
    output:
        aims = f"results/ref/all_AIMs_{genome1}_{genome2}",
        current_aims_file = "results/ref/current_aims_file",
        formatted = f"{aims}_formatted",
        aims_bed = f"{aims}.mod",
        aims_truebed = f"{aims}.mod.bed",
    shell:
        "cp {input.provided_aims} {output.aims} ; "
        "echo {output.aims} > {output.current_aims_file} ; "
        "perl -pi -e 's/_/-/g' {output.aims} ; "
        "cp {output.aims} {output.formatted}; "
        "perl -pi -e 's/\t/_/g' {output.aims} ; "
        """awk -v OFS='\t' '{{print $1"_"$2,$2,$3,$4}}' {output.formatted} > {output.aims_bed} ; """  
        """ awk -v OFS='\t' '{{print $1,$2,$2,$3,$4}}' {output.formatted} > {output.aims_truebed}  """
