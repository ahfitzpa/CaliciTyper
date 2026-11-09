#!/bin/bash -l

# Set the number of nodes

#SBATCH -N 1

# Set the number of tasks/cores per node required 
#SBATCH -n 3

# Set the walltime of the job to 1 hour (format is hh:mm:ss)
#SBATCH -t 01:00:00

# E-mail on begin (b), abort (a) and end (e) of job
#SBATCH --mail-type=ALL

# E-mail address of recipient
#SBATCH --mail-user=myemailaddress@ucd.ie

# Specifies the jobname
#SBATCH --job-name=myjob

# Change working directory to current directory
# liuyang open if submitted to batch serve !!!
# cd $SLURM_SUBMIT_DIR
# Your code here! (This example just prints the current date and exits!)

python "$(dirname "$0")/makedb.py" \
       "$(dirname "$0")/../norotyper/data/database/refseq_2.fasta" \
       "$(dirname "$0")/../norotyper/data/database/refseqdb_2"

