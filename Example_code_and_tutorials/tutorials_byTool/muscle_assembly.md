# Muscle

```

#!/bin/bash
#SBATCH -J muscle_run -o muscle.out -t 20:00:00 --nodes=1 --ntasks-per-node=1 --cpus-per-task=12 --mem=50G
# Remember to change the "queueSize" parameter in the config/*_slurm.config file that you are using. This corresponds with "--ntasks" 
# to control how many jobs are submitted at once.
# This script works on TAMU's HPRC, but you need to follow the instructions on the Github to get the right conda 
# environment installed on your computing environment

module load GCCcore/11.2.0 MUSCLE/5.1

muscle -super5 concat_RespSwab_290b_PCV_9_extracted.fasta -output concat_muscle_output

```