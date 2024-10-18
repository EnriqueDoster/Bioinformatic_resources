#!/bin/bash
#SBATCH -J dada2_3 -t 01:00:00 -n 1 --cpus-per-task=12 -o 3_dada2.out

# initiate R
R --vanilla

# run Rscript (3_dada2.R)
Rscript 3_dada2.R
