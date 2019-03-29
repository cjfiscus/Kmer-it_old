#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --output=%j.stdout
#SBATCH --error=%j.stderr
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH --job-name="test"
#SBATCH -p koeniglab

axel ftp.sra.ebi.ac.uk/vol1/fastq/SRR194/007/SRR1945447/SRR1945447.fastq.gz -o ./temp/470.1/470.1.fastq.gz  
