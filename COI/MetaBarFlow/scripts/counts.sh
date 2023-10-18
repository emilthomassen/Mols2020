#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 8G
#SBATCH -c 2
#SBATCH --time 2:00:00

Rscript counts_coleoptera.r
