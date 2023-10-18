#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 8G
#SBATCH -c 2
#SBATCH --time 2:00:00

Rscript merge_files_cluster_curated_data.r
