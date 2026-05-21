#!/bin/bash
#
# One-shot setup for the GPU Training on CoreHPC workshop.
# Run this on the CoreHPC login node from any directory.
#
# What it does:
#   1. Picks a workshop root under /mnt/scratch/user/$USER
#   2. Clones karpathy/nanoGPT into that root
#   3. Creates a conda env from environment.yml
#   4. Runs the Shakespeare char-level data prep step
#   5. Tells you which sbatch script to run first

set -euo pipefail

if [ -z "${USER:-}" ]; then
    echo "ERROR: \$USER is not set."
    exit 1
fi

WORKSHOP_ROOT="/mnt/scratch/user/${USER}/gpu-training-on-corehpc"
NANOGPT_DIR="${WORKSHOP_ROOT}/nanoGPT"
ENV_NAME="gpu-training-corehpc"
ENV_FILE="$(cd "$(dirname "$0")" && pwd)/environment.yml"

echo "Workshop root: ${WORKSHOP_ROOT}"
mkdir -p "${WORKSHOP_ROOT}"

# 1. Clone nanoGPT (fresh)
if [ -d "${NANOGPT_DIR}/.git" ]; then
    echo "nanoGPT already cloned at ${NANOGPT_DIR}, pulling latest."
    git -C "${NANOGPT_DIR}" pull --ff-only
else
    echo "Cloning nanoGPT into ${NANOGPT_DIR}"
    git clone https://github.com/karpathy/nanoGPT.git "${NANOGPT_DIR}"
fi

# 2. Load miniforge and build the conda env
if ! command -v module &>/dev/null; then
    echo "ERROR: 'module' command not found. Are you on a CoreHPC node?"
    exit 1
fi

echo "Loading miniforge3 module."
module load miniforge3

# Redirect conda envs/pkgs onto scratch (home is only 20 GB on CoreHPC)
conda config --append envs_dirs  "/mnt/scratch/user/${USER}/conda/envs" 2>/dev/null || true
conda config --append pkgs_dirs  "/mnt/scratch/user/${USER}/conda/pkgs" 2>/dev/null || true

if conda env list | awk '{print $1}' | grep -qx "${ENV_NAME}"; then
    echo "Conda env '${ENV_NAME}' already exists. Skipping create."
else
    echo "Creating conda env '${ENV_NAME}' from ${ENV_FILE}"
    conda env create -n "${ENV_NAME}" -f "${ENV_FILE}"
fi

# 3. Prepare the Shakespeare char-level dataset (downloads ~1 MB of text)
echo "Preparing the Shakespeare char-level dataset."
# shellcheck source=/dev/null
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${ENV_NAME}"
python "${NANOGPT_DIR}/data/shakespeare_char/prepare.py"

# 4. Print what to do next
cat <<EOF

Setup complete.

  nanoGPT clone:    ${NANOGPT_DIR}
  conda env:        ${ENV_NAME}
  dataset:          ${NANOGPT_DIR}/data/shakespeare_char/{train,val}.bin

To run the first job, from this directory:

  sbatch slurm/01_single_gpu.sbatch

Then walk through notes/lab_guide.md.

EOF
