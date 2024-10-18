#!/bin/bash
#SBATCH -J dada2_2 -t 01:00:00 -n 1 --cpus-per-task=12 -o 2_dada2.out

# initiate R
R --vanilla

# run Rscript (2_dada2.R)
Rscript 2_dada2.R
