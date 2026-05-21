#!/bin/bash
#
# 03_multi_gpu_ddp.sbatch
# ------------------------
# DDP across 2 GPUs on a single node using torchrun. Requires the `pod`
# partition because small_gpu / large_gpu cap at 1 GPU per running job.
#
# Submit ONLY if your CoreHPC account is approved for pod access:
#   sbatch slurm/03_multi_gpu_ddp.sbatch
#
# Otherwise read this script as a reference for the structure of a multi-GPU
# job; the SBATCH lines and the torchrun invocation are the parts that matter.
#
#SBATCH --job-name=gpu-train-03-ddp
#SBATCH --partition=pod
#SBATCH --gres=gpu:nvidia_h200:2
#SBATCH --cpus-per-task=8
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH --time=00:30:00
#SBATCH --output=logs/03-%j.out
#SBATCH --error=logs/03-%j.err

set -euo pipefail
mkdir -p logs

WORKSHOP_ROOT="/mnt/scratch/user/${USER}/gpu-training-on-corehpc"
NANOGPT_DIR="${WORKSHOP_ROOT}/nanoGPT"
OUT_DIR="${WORKSHOP_ROOT}/runs/${SLURM_JOB_ID}"
mkdir -p "${OUT_DIR}"

module load miniforge3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate gpu-training-corehpc

cd "${NANOGPT_DIR}"

# --nproc_per_node MUST equal the GPU count from --gres above.
# --standalone is single-node only; for multi-node use --master_addr / --master_port.
torchrun --standalone --nproc_per_node=2 \
    train.py config/train_shakespeare_char.py \
    --out_dir="${OUT_DIR}" \
    --compile=False

echo
echo "Multi-GPU training complete."
echo "Checkpoint: ${OUT_DIR}/ckpt.pt"
