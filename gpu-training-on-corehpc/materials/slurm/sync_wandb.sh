#!/bin/bash
#
# sync_wandb.sh
# --------------
# Run this on the CoreHPC LOGIN node (not a compute node) after an offline
# W&B run has completed. It will upload all unsynced offline-run-* directories
# to your W&B dashboard.
#
# First time only: run `wandb login` to paste your API key.
#
# Usage:
#   bash slurm/sync_wandb.sh
#

set -euo pipefail

WORKSHOP_ROOT="/mnt/scratch/user/${USER}/gpu-training-on-corehpc"
WANDB_DIR="${WORKSHOP_ROOT}/wandb"

if [ ! -d "${WANDB_DIR}" ]; then
    echo "No W&B directory found at ${WANDB_DIR}. Did the offline job run?"
    exit 1
fi

module load miniforge3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate gpu-training-corehpc

# Sync all unsynced offline runs
wandb sync "${WANDB_DIR}"/offline-run-*
