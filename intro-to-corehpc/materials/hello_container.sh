#!/bin/bash
#SBATCH --job-name=hello_container
#SBATCH --partition=cpu
#SBATCH --time=00:05:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --output=hello_container_%j.out
#SBATCH --error=hello_container_%j.err

SIF=hello-world_1.0.sif
WORKDIR=/mnt/scratch/user/$USER/container_demo

mkdir -p "$WORKDIR"

# --bind <path on the cluster>:<path inside the container>
srun apptainer exec --bind "$WORKDIR":/data "$SIF" \
    bash -c 'hi > /data/greeting.txt'

echo "--- wrote $WORKDIR/greeting.txt:"
cat "$WORKDIR/greeting.txt"
