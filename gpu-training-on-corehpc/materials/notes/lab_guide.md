# Lab guide

Step-by-step walk-through of the four sbatch examples. The video covers the *why* for each pattern; this guide is the *what to type*.

Before starting, make sure you've run through `pre_workshop_setup.md` and that `bash setup.sh` has finished.

You should be on the **CoreHPC login node** (`chpc-gs-login-vm1`) in the `materials/` directory.

---

## Exercise 1: baseline single-GPU run

**Goal:** confirm your environment works and you know how to inspect a running job.

```bash
sbatch slurm/01_single_gpu.sbatch
squeue -u $USER
```

Note the JobID from `squeue`. While the job is running:

```bash
# One-shot snapshot of GPU utilization on the job's node:
srun --jobid=<JOBID> --pty nvidia-smi
```

You should see GPU-Util hovering near 100% and Memory-Usage growing.

When the job finishes (~3-5 minutes):

```bash
seff <JOBID>
sacct -j <JOBID> --format=JobID,State,Elapsed,MaxRSS,ReqTRES,AllocTRES

# Inspect what got saved:
ls -lh /mnt/scratch/user/$USER/gpu-training-on-corehpc/runs/<JOBID>/
# Should contain: ckpt.pt, gpu_util.log
```

**What to look for:** `ReqTRES` and `AllocTRES` agree that you got 1 L40s GPU. The first column of `gpu_util.log` is timestamps; the second is GPU-Util percentage over time.

---

## Exercise 2: same run with offline W&B logging

**Goal:** see how W&B works when the compute node has no internet.

```bash
sbatch slurm/02_wandb_offline.sbatch
squeue -u $USER
```

When the job finishes:

```bash
ls /mnt/scratch/user/$USER/gpu-training-on-corehpc/wandb/
# Should contain: offline-run-YYYYMMDD_HHMMSS-<id>/

# Sync to your W&B dashboard (login node has internet via proxy):
bash slurm/sync_wandb.sh
```

Open your project at <https://wandb.ai/$your_username/gpu-training-corehpc>. You should see:

-   The training/validation loss curves
-   **System metrics** tab: per-device GPU utilization, memory, power, temperature
-   The run name should be `run-<SLURM_JOB_ID>`

**What to look for in the dashboard:** GPU-Util on device 0 should sit near 100% during training. There's only one device because you asked for `--gres=...:1`. This is the visualization that makes over-reservation obvious; see Exercise 3 for the contrast.

---

## Exercise 3: multi-GPU DDP (pod partition only)

**Goal:** see how `torchrun` distributes a single training run across multiple GPUs.

If you have `pod` access:

```bash
sbatch slurm/03_multi_gpu_ddp.sbatch
squeue -u $USER
```

While running, check both GPUs:

```bash
srun --jobid=<JOBID> --pty nvidia-smi
```

Both devices should show high utilization. Both processes contribute to the same effective batch.

If you don't have `pod` access, **read** `slurm/03_multi_gpu_ddp.sbatch` and notice:

-   `--gres=gpu:nvidia_h200:2` reserves 2 GPUs
-   `torchrun --nproc_per_node=2` launches 2 processes, one per GPU
-   The two numbers **must match**. If you set `--nproc_per_node=1` but `--gres=...:2`, you're paying for 2 GPUs and only using 1. This is the over-reservation trap from the video.

**Optional follow-up:** modify the script to set `--nproc_per_node=1` while keeping `--gres=...:2`, submit it, and look at the W&B system metrics tab. You'll see device 0 at 100% and device 1 at 0% the whole run.

---

## Exercise 4: graceful preemption + resume

**Goal:** prove that the checkpointing + requeue setup actually recovers from interruption.

First run (starts from scratch, writes a checkpoint, runs to completion):

```bash
sbatch slurm/04_preempt_requeue.sbatch
squeue -u $USER
```

Watch for the message `[..] no checkpoint, starting from scratch.` in `logs/04-<JOBID>.out`.

After the job finishes, look at the checkpoint:

```bash
ls -lh /mnt/scratch/user/$USER/gpu-training-on-corehpc/runs/persistent/
# Should contain ckpt.pt
```

Now **simulate a preemption** by killing a fresh submission with SIGTERM mid-run, then resubmitting:

```bash
# Resubmit:
sbatch slurm/04_preempt_requeue.sbatch
squeue -u $USER

# Wait ~30 seconds for it to actually be running, then:
scancel --signal=TERM <JOBID>

# Look at the log. You should see the SIGTERM trap fire:
tail logs/04-<JOBID>.out
# [<timestamp>] received SIGTERM, exiting cleanly. ckpt is at .../persistent/ckpt.pt

# Now resubmit. The new job should resume from the checkpoint:
sbatch slurm/04_preempt_requeue.sbatch
tail -f logs/04-<NEW_JOBID>.out
# [<timestamp>] found existing checkpoint, resuming.
```

**What to look for:** the loss in the resumed run should pick up roughly where the previous run left off, not from the initial value. The optimizer state is restored along with the model state, so there's no spike from a fresh Adam initialization.

---

## Cleanup

Workshop materials live under `/mnt/scratch/user/$USER/gpu-training-on-corehpc/`. Scratch is 30-day auto-cleanup (after GA) and will be wiped before GA. If you want to keep any of the trained checkpoints, copy them to your lab's HIVE share:

```bash
cp /mnt/scratch/user/$USER/gpu-training-on-corehpc/runs/persistent/ckpt.pt \
   /mnt/gladstone/<lab>/models/shakespeare-char-ckpt.pt
```

HIVE (`/mnt/gladstone`) is the only backed-up storage on CoreHPC.
