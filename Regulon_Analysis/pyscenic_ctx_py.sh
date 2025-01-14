#!/bin/bash

#SBATCH --job-name="pyscenic_ctx"       # <--- Assigning a job is not required, but recommended; it is for your reference
#SBATCH --mem=200g                    # <--- This job will require 20 cores
#SBATCH -A bms_q                    # <--- Run in the genacc_q Slurm account; this is the default account if none specified
#SBATCH --mail-type="ALL"       # <--- Email your official FSU email address on all job events


source /gpfs/home/lg23w/anaconda3/etc/profile.d/conda.sh
conda activate pyscenic

pyscenic ctx \
    /gpfs/home/lg23w/lipoglucotoxicity/SCENIC/lipo_integrated_endocrine.adj.tsv \
    /gpfs/home/lg23w/lipoglucotoxicity/SCENIC/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
    /gpfs/home/lg23w/lipoglucotoxicity/SCENIC/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
    --annotations_fname /gpfs/home/lg23w/lipoglucotoxicity/SCENIC/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl \
    --expression_mtx_fname /gpfs/home/lg23w/lipoglucotoxicity/SCENIC/lipo_integrated_endocrine.loom \
    --mode "dask_multiprocessing" \
    --output /gpfs/home/lg23w/lipoglucotoxicity/SCENIC/lipo_integrated_endocrine_regulons.csv \
    --num_workers 51