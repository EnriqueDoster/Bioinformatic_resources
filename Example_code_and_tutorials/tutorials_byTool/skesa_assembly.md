# Skesa

Works on grace to create contigs.

```


#!/bin/bash
#SBATCH -J skesa_run -o skesa.out -t 20:00:00 --nodes=1 --ntasks-per-node=1 --cpus-per-task=12 --mem=50G
# Remember to change the "queueSize" parameter in the config/*_slurm.config file that you are using. This corresponds with "--ntasks" 
# to control how many jobs are submitted at once.
# This script works on TAMU's HPRC, but you need to follow the instructions on the Github to get the right conda 
# environment installed on your computing environment

module load GCC/7.3.0-2.30  OpenMPI/3.1.1 SKESA/2.3.0

skesa --fasta concat_RespSwab_290b_PCV_9_extracted.fasta --contigs_out Skesa_PCV_9_contigs --cores 12

```