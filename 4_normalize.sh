#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=400gb
#SBATCH --output=%j.stdout
#SBATCH --error=%j.stderr
#SBATCH --time=2-00:00:00
#SBATCH --mail-user=cfisc004@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="norm"
#SBATCH -p koeniglab

# count tables 
INPUT=$1
GCNORM="temp.txt"
GCTMMNORM=$2
PROP_TABLE="prop.txt"

# gc correct counts
Rscript 4a_gc_correct.R "$INPUT" "$GCNORM" "$PROP_TABLE"  

# tmm normalize counts
Rscript 4b_tmm.R "$GCNORM" "$GCTMMNORM"
