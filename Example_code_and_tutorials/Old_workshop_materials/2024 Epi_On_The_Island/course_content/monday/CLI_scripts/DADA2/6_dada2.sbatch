#!/bin/bash
#SBATCH -J dada2_6 -t 01:00:00 -n 1 --cpus-per-task=1 -o 6_dada2.out

# initiate R
R --vanilla

# run Rscript (6_dada2.R)
Rscript 6_dada2.R
