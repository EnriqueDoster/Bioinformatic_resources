#!/bin/bash
#SBATCH -J dada2_7 -t 01:00:00 -n 1 --cpus-per-task=12 -o 7_dada2.out

# initiate R
R --vanilla

# run Rscript (7_dada2.R)
Rscript 7_dada2.R
