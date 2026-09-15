#!/bin/bash
#SBATCH -p medmem,standard 
#SBATCH --job-name=ancestryhmm
#SBATCH -t 24:00:00
#SBATCH -c 1 

ancestry_hmm -a 2 0.5 0.5 -p 0 -10000 0.5 -p 1 30 0.5 -e 0.05 -s current.samples.list_allchrs -i all.indivs.hmm.combined_allchrs.filtered_focalchroms_NC-048565.1
#ancestry_hmm -a 2 0.5 0.5 -p 0 -10000 0.5 -p 1 30 0.5 -e 0.05 -s current.samples.list_allchrs -i all.indivs.hmm.combined_allchrs.filtered_focalchroms
