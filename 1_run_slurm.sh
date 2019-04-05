#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=16G
#SBATCH --output=%j.stdout
#SBATCH --error=%j.stderr
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH --job-name="test"
#SBATCH -p koeniglab
#SBATCH --array=2-6

# software versions
#samtools 1.8; trimmomatic 0.36; bedtools 2.27.0; jellyfish 2.2.9; bwa 0.7.17

# load required modules (slurm)
module load trimmomatic/0.36 jellyfish/2.2.9 samtools/1.9 bwa/0.7.17

# run software
sh kmerit.sh params "$SLURM_ARRAY_TASK_ID"
