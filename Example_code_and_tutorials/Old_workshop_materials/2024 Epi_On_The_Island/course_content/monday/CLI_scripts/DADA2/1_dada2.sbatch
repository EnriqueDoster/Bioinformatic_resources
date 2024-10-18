#!/bin/bash
#SBATCH -J dada2_1 -t 01:00:00 -n 1 --cpus-per-task=1 -o 1_dada2.out

# initiate R
R --vanilla

# run Rscript (1_dada2.R)
Rscript 1_dada2.R
