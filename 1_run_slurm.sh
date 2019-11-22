#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=16G
#SBATCH --output=std/%j.stdout
#SBATCH --error=std/%j.stderr
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="sur"
#SBATCH -p koeniglab
#SBATCH --array=2-20

# software versions
#samtools 1.8; trimmomatic 0.36; bedtools 2.27.0; jellyfish 2.2.9; bwa 0.7.17

# load required modules (slurm)
module load trimmomatic/0.36 jellyfish/2.2.9 samtools/1.9 bwa/0.7.17 

# run software
sh scripts/kmerit.sh ./params3 "$SLURM_ARRAY_TASK_ID"
