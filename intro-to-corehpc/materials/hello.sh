#!/bin/bash
#
# Minimal CoreHPC example job. Submit it from a login node with:
#   sbatch hello.sh
# Output lands in hello_<jobid>.out in the directory you submit from.
#
#SBATCH --job-name=hello
#SBATCH --partition=cpu
#SBATCH --time=00:05:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --output=hello_%j.out   # %j = job ID
#SBATCH --error=hello_%j.err

srun echo "Hello from $(hostname)"
