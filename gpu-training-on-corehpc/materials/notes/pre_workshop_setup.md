# Pre-workshop setup

Do these before pressing play on the video so you can run the hands-on parts in real time.

## 1. Confirm your CoreHPC account works

You should be able to log in to the Gladstone-side CoreHPC cluster:

```bash
# From the Gladstone network (wired, Gladstone-MB, or Ivanti VPN), on a Gladstone device:
ssh <user>@10.98.160.34          # chpc-gs-bastion-vm1 (address is temporary)
# Enter your UCSF MyAccess password, approve the DUO push.
ssh chpc-gs-login-vm1
```

If you don't have an account yet, request one via the [CoreHPC support form](https://tiny.ucsf.edu/Gethpc) and wait for confirmation before doing the rest of these steps. See the [Wynton to CoreHPC migration guide](../../../intro-to-corehpc/Wynton_to_CoreHPC_Migration.qmd) for the full bastion + login walk-through.

## 2. Request pod access (optional)

The multi-GPU example (`03_multi_gpu_ddp.sbatch`) requires the `pod` partition, which is request-only. If you want to run it during the workshop:

-   Submit a request via the [CoreHPC support form](https://tiny.ucsf.edu/Gethpc)
-   Mention what you want to train and roughly how much GPU time you need

If you don't have pod access, you can still read the script and submit `01`, `02`, and `04` (all on `small_gpu`).

## 3. Set up a Weights & Biases account

-   Sign up at [wandb.ai](https://wandb.ai)
-   Once logged in, grab your API key from <https://wandb.ai/authorize>
-   On the CoreHPC login node, after `setup.sh` has built the conda env:

```bash
module load miniforge3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate gpu-training-corehpc
wandb login    # paste your API key
```

`wandb login` writes a credentials file to your home directory. You only need to do it once per account.

## 4. (Optional) Install Globus Connect Personal

If you want to pull trained model weights off CoreHPC at the end of the workshop, set up [Globus Connect Personal](https://www.globus.org/globus-connect-personal) on your laptop. Lab Hive shares (`/mnt/gladstone/<lab>/`) are the recommended long-term home for anything you can't afford to lose. Scratch is volatile.

## 5. Run the setup script

```bash
# On the CoreHPC login node:
cd /mnt/scratch/user/$USER
git clone https://github.com/gladstone-institutes/Bioinformatics-Workshops.git
cd Bioinformatics-Workshops/gpu-training-on-corehpc/materials
bash setup.sh
```

Takes ~5 minutes the first time (mostly conda package downloads). On rerun it's a no-op.

When `setup.sh` finishes you're ready to follow `lab_guide.md`.
