#!/bin/bash
#SBATCH -p medmem,standard 
#SBATCH --job-name=perl_launch
#SBATCH -t 2:00:00
#SBATCH -c 8 

perl debug.pl
