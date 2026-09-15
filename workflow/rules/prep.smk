
# rename chromosome names to replace underscores with dashes
# unfortunately cannot have underscores in chrom names
rule rename_chrom_names:
    input: 
        og_g1 = og_genome1,
        og_g2 = og_genome2,
    output: 
        renamed_g1 = og_genome1+"_chr-renamed",
        renamed_g2 = og_genome2+"_chr-renamed",
    shell:
        "perl -p -e 's/_/-/g' {input.og_g1} > {output.renamed_g1} ; "
        "perl -p -e 's/_/-/g' {input.og_g2} > {output.renamed_g2} "

# reformat the AIMs file for downstream compatibility
# make "current_aims_file"
rule reformat_aims:
    input:
        provided_aims=config["provide_AIMs"],
    output:
        aims = "all_AIMs_"+genome1+"_"+genome2,
        current_aims_file="current_aims_file",
        formatted = aims+"_formatted",
        aims_bed = aims+".mod",
        aims_truebed=aims+".mod.bed",
    shell:
        "cp {input.provided_aims} {output.aims} ; "
        "echo {output.aims} > {output.current_aims_file} ; "
        "perl -pi -e 's/_/-/g' {output.aims} ; "
        "cp {output.aims} {output.formatted}; "
        "perl -pi -e 's/\t/_/g' {output.aims} ; "
        """awk -v OFS='\t' '{{print $1"_"$2,$2,$3,$4}}' {output.formatted} > {output.aims_bed} ; """  
        """ awk -v OFS='\t' '{{print $1,$2,$2,$3,$4}}' {output.formatted} > {output.aims_truebed}  """
