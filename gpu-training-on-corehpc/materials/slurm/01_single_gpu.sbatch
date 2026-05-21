#!/bin/bash
#
# 01_single_gpu.sbatch
# ---------------------
# Baseline: one L40s, no W&B, no DDP. Trains the Shakespeare char-level model
# in ~3-5 minutes. Demonstrates the minimum-viable SLURM script and the
# checkpoint file landing on scratch.
#
# Submit with:
#   sbatch slurm/01_single_gpu.sbatch
#
# Watch with:
#   squeue -u $USER
#   srun --jobid=<JOBID> --pty nvidia-smi    # one-shot util snapshot
#
#SBATCH --job-name=gpu-train-01
#SBATCH --partition=small_gpu
#SBATCH --gres=gpu:nvidia_l40s:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --output=logs/01-%j.out
#SBATCH --error=logs/01-%j.err

set -euo pipefail
mkdir -p logs

# Paths set up by setup.sh
WORKSHOP_ROOT="/mnt/scratch/user/${USER}/gpu-training-on-corehpc"
NANOGPT_DIR="${WORKSHOP_ROOT}/nanoGPT"
OUT_DIR="${WORKSHOP_ROOT}/runs/${SLURM_JOB_ID}"
mkdir -p "${OUT_DIR}"

# Activate the conda env created by setup.sh
module load miniforge3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate gpu-training-corehpc

# Background nvidia-smi loop so we have a record of GPU util after the job ends
nvidia-smi --query-gpu=timestamp,utilization.gpu,memory.used,memory.total \
           --format=csv --loop=5 > "${OUT_DIR}/gpu_util.log" &
NVSMI_PID=$!
trap 'kill ${NVSMI_PID} 2>/dev/null || true' EXIT

cd "${NANOGPT_DIR}"
python train.py config/train_shakespeare_char.py \
    --out_dir="${OUT_DIR}" \
    --compile=False

echo
echo "Training complete."
echo "Checkpoint: ${OUT_DIR}/ckpt.pt"
echo "GPU log:    ${OUT_DIR}/gpu_util.log"
echo
echo "After this job finishes, run:"
echo "  seff ${SLURM_JOB_ID}"
echo "  sacct -j ${SLURM_JOB_ID} --format=JobID,State,Elapsed,MaxRSS,ReqTRES,AllocTRES"
