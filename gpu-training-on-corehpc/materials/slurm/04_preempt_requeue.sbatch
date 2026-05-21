#!/bin/bash
#
# 04_preempt_requeue.sbatch
# --------------------------
# Demonstrates graceful preemption: --requeue tells SLURM to resubmit the job
# if it gets preempted, the SIGTERM trap gives us a clean shutdown window,
# and init_from=resume makes the second run pick up from the last checkpoint.
#
# Submit with:
#   sbatch slurm/04_preempt_requeue.sbatch
#
# To simulate a preemption (since we can't trigger a real one on demand):
#   scancel --signal=TERM <JOBID>     # send SIGTERM
#   sbatch slurm/04_preempt_requeue.sbatch  # resubmit, will resume from ckpt
#
# Note: the OUT_DIR below is fixed (not job-id scoped) so the resumed run
# can find the previous checkpoint.
#
#SBATCH --job-name=gpu-train-04-requeue
#SBATCH --partition=small_gpu
#SBATCH --gres=gpu:nvidia_l40s:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --requeue
#SBATCH --signal=B:SIGTERM@60
#SBATCH --output=logs/04-%j.out
#SBATCH --error=logs/04-%j.err

set -euo pipefail
mkdir -p logs

WORKSHOP_ROOT="/mnt/scratch/user/${USER}/gpu-training-on-corehpc"
NANOGPT_DIR="${WORKSHOP_ROOT}/nanoGPT"
OUT_DIR="${WORKSHOP_ROOT}/runs/persistent"      # fixed dir so resume can find the checkpoint
mkdir -p "${OUT_DIR}"

# Trap SIGTERM: print a marker, let the script exit cleanly.
# The checkpoint is already on disk from the last eval interval.
trap 'echo "[$(date)] received SIGTERM, exiting cleanly. ckpt is at ${OUT_DIR}/ckpt.pt"; exit 0' SIGTERM

module load miniforge3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate gpu-training-corehpc

# Decide init_from based on whether a checkpoint already exists
if [ -f "${OUT_DIR}/ckpt.pt" ]; then
    INIT_FROM="resume"
    echo "[$(date)] found existing checkpoint, resuming."
else
    INIT_FROM="scratch"
    echo "[$(date)] no checkpoint, starting from scratch."
fi

cd "${NANOGPT_DIR}"
python train.py config/train_shakespeare_char.py \
    --out_dir="${OUT_DIR}" \
    --init_from="${INIT_FROM}" \
    --always_save_checkpoint=True \
    --compile=False &

# `wait` is what lets the SIGTERM trap fire while python is running.
wait $!
