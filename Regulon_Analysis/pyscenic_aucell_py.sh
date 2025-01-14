#!/bin/bash

#SBATCH --job-name="pyscenic_aucell"       # <--- Assigning a job is not required, but recommended; it is for your reference
#SBATCH --mem=200g                    # <--- This job will require 20 cores
#SBATCH -A bms_q                    # <--- Run in the genacc_q Slurm account; this is the default account if none specified
#SBATCH --mail-type="ALL"       # <--- Email your official FSU email address on all job events


source /gpfs/home/lg23w/anaconda3/etc/profile.d/conda.sh
conda activate pyscenic
 

pyscenic aucell \
    --num_workers 51 \
    -o /gpfs/home/lg23w/lipoglucotoxicity/SCENIC/lipo_integrated_endocrine_SCENIC.loom \
    /gpfs/home/lg23w/lipoglucotoxicity/SCENIC/lipo_integrated_endocrine_new.loom \
    /gpfs/home/lg23w/lipoglucotoxicity/SCENIC/lipo_integrated_endocrine_regulons.csv