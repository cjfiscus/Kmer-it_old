#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=16G
#SBATCH --output=%j.stdout
#SBATCH --error=%j.stderr
#SBATCH --job-name="kmerit"
#SBATCH -p batch
#SBATCH --array=2-20

# software versions
#samtools 1.9; trimmomatic 0.36; bedtools 2.29.2; jellyfish 2.2.9; bwa 0.7.17

# load required modules (slurm)
module load trimmomatic/0.36 jellyfish/2.2.9 samtools/1.9 bwa/0.7.17 

# run software
## $SLURM_ARRAY_TASK_ID is the sample to run from sample.txt
sh scripts/kmerit.sh ./params "$SLURM_ARRAY_TASK_ID"
