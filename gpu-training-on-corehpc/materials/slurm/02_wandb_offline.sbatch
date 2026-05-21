#!/bin/bash
#
# 02_wandb_offline.sbatch
# ------------------------
# Same training run as 01, but with Weights & Biases logging in offline mode.
# Compute nodes on CoreHPC have no internet, so we log to disk during the run
# and sync from the login node afterwards.
#
# Submit with:
#   sbatch slurm/02_wandb_offline.sbatch
#
# After the job finishes, from the LOGIN node:
#   bash slurm/sync_wandb.sh
#
#SBATCH --job-name=gpu-train-02
#SBATCH --partition=small_gpu
#SBATCH --gres=gpu:nvidia_l40s:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --output=logs/02-%j.out
#SBATCH --error=logs/02-%j.err

set -euo pipefail
mkdir -p logs

WORKSHOP_ROOT="/mnt/scratch/user/${USER}/gpu-training-on-corehpc"
NANOGPT_DIR="${WORKSHOP_ROOT}/nanoGPT"
OUT_DIR="${WORKSHOP_ROOT}/runs/${SLURM_JOB_ID}"
mkdir -p "${OUT_DIR}"

module load miniforge3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate gpu-training-corehpc

# W&B in offline mode. Sync happens later from the login node.
export WANDB_MODE=offline
export WANDB_DIR="${WORKSHOP_ROOT}/wandb"
export WANDB_NAME="run-${SLURM_JOB_ID}"
export WANDB_PROJECT="gpu-training-corehpc"
mkdir -p "${WANDB_DIR}"

cd "${NANOGPT_DIR}"
python train.py config/train_shakespeare_char.py \
    --out_dir="${OUT_DIR}" \
    --compile=False \
    --wandb_log=True \
    --wandb_project="${WANDB_PROJECT}" \
    --wandb_run_name="${WANDB_NAME}"

echo
echo "Training complete."
echo "W&B offline run saved under: ${WANDB_DIR}"
echo
echo "From the LOGIN node, sync to W&B:"
echo "  bash $(dirname "$0")/sync_wandb.sh"
