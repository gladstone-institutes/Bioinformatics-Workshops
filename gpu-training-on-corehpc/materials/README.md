# GPU Training on CoreHPC: hands-on materials

Companion materials for the **GPU Training on CoreHPC** workshop video. The example workload is Karpathy's [nanoGPT](https://github.com/karpathy/nanoGPT) (MIT license) training a character-level Shakespeare model. It's small enough to finish in a few minutes on one L40s and ships with built-in checkpointing, DDP, and W&B integration, so all four pitfalls covered in the video can be demonstrated on a real run.

## What's in here

```
materials/
├── README.md              this file
├── setup.sh               one-shot: clones nanoGPT, builds the conda env, prepares the dataset
├── environment.yml        conda environment (torch, tiktoken, wandb, numpy)
├── slurm/
│   ├── 01_single_gpu.sbatch        baseline run on one L40s
│   ├── 02_wandb_offline.sbatch     same run + offline W&B logging
│   ├── 03_multi_gpu_ddp.sbatch     torchrun across 2 H200s (pod partition)
│   ├── 04_preempt_requeue.sbatch   --requeue + SIGTERM handler + resume
│   └── sync_wandb.sh               one-liner to sync offline runs from the login node
└── notes/
    ├── pre_workshop_setup.md       what to do before playing the video
    └── lab_guide.md                step-by-step walkthrough
```

## Prerequisites

-   CoreHPC account ([request via UCSF](https://tiny.ucsf.edu/Gethpc))
-   On the Gladstone network (wired, Gladstone-MB, or Ivanti VPN) on a Gladstone device
-   A W&B account ([wandb.ai](https://wandb.ai))
-   `pod` partition access if you want to run script `03`; ask via the support form

## Quick start

```bash
# Log in to CoreHPC (see the migration deck for the bastion hop)
ssh <user>@10.98.160.34
ssh chpc-gs-login-vm1

# Clone this repo to scratch and run setup
cd /mnt/scratch/user/$USER
git clone https://github.com/gladstone-institutes/Bioinformatics-Workshops.git
cd Bioinformatics-Workshops/gpu-training-on-corehpc/materials
bash setup.sh

# Submit your first job
sbatch slurm/01_single_gpu.sbatch
squeue -u $USER
```

Then follow `notes/lab_guide.md`.
