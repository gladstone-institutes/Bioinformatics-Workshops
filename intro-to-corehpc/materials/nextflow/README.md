# Nextflow on CoreHPC: templates

Two project-agnostic files for running Nextflow pipelines on CoreHPC's SLURM
scheduler. Copy them into your project and set one environment variable; you
should not need to edit either one.

| File | Purpose |
|---|---|
| `corehpc_nextflow_template.config` | SLURM executor, resource ceilings, containers, reports |
| `launch_corehpc_template.sh` | Submits Nextflow as a driver job; has a built-in smoke test |
| `running_nextflow_on_corehpc.md` | The original write-up the full guide was adapted from |

The login node kills your processes when the SSH session drops, so Nextflow
itself is submitted as a batch job, "the driver". It sits on a compute node and
submits one child job per task.

## Quick start

```bash
export COREHPC_PROJECT=/mnt/gladstone/<pi>/<share>/<project>

./launch_corehpc_template.sh smoke   # end-to-end check, no real data
./launch_corehpc_template.sh check   # print the driver sbatch, submit nothing
```

A real run needs a local copy of the pipeline and a samplesheet:

```bash
PIPELINE=$COREHPC_PROJECT/assets/nf-core-rnaseq-3.19.0/3_19_0 \
SAMPLESHEET=$COREHPC_PROJECT/input/samplesheet.csv \
  ./launch_corehpc_template.sh
```

Both files take optional overrides from the environment: `COREHPC_ACCOUNT`,
`COREHPC_PARTITION`, `COREHPC_WORKDIR`, `COREHPC_CACHEDIR`, `DRIVER_TIME`,
`DRIVER_MEM`.

## Before the first run

Compute nodes have no internet. On the **login node** you have to install
Nextflow and Java, then pre-fetch three things: the pipeline code, its container
images, and the plugins it pins. The guide walks through each.

## Full guide

[Running Nextflow on CoreHPC](https://gladstone-institutes.github.io/Bioinformatics-Workshops/Intro_to_CoreHPC/Nextflow_on_CoreHPC.html)
covers installation, the pre-fetches, surviving the driver walltime, a
pre-flight checklist, and an error quick reference.

Scratch is `/mnt/scratch` on both halves. Persistent storage differs: these
examples use Gladstone's `/mnt/gladstone`, so substitute your FAC share on the
UCSF half.
