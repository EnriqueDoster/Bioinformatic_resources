#!/bin/bash
#SBATCH -J dada2_4 -t 01:00:00 -n 1 --cpus-per-task=12 -o 4_dada2.out

# initiate R
R --vanilla

# run Rscript (4_dada2.R)
Rscript 4_dada2.R
