#!/bin/bash

#SBATCH --job-name="pyscenic_grn"       # <--- Assigning a job is not required, but recommended; it is for your reference
#SBATCH --mem=250g                    # <--- This job will require 20 cores
#SBATCH -A bms_q                    # <--- Run in the genacc_q Slurm account; this is the default account if none specified
#SBATCH --mail-type="ALL"       # <--- Email your official FSU email address on all job events


source /gpfs/home/lg23w/anaconda3/etc/profile.d/conda.sh
conda activate pyscenic
#module load R/4.3.2      # <--- Load the 'R/4.3.2' module\
 


pyscenic grn \
    --num_workers 63 \
    -o /gpfs/home/lg23w/lipoglucotoxicity/SCENIC/lipo_integrated.adj.tsv \
    /gpfs/home/lg23w/lipoglucotoxicity/SCENIC/lipo_integrated.loom\
    /gpfs/home/lg23w/lipoglucotoxicity/SCENIC/allTFs_hg38.txt
