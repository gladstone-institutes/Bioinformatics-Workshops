# Running Nextflow / nf-core pipelines on UCSF CoreHPC

A setup guide for running nf-core pipelines on CoreHPC. Written from the
GB-TH-1359 step6 run (nf-core/rnaseq, 56 PBMC samples); everything here was
learned by hitting the failure first.

**Start here.** Two project-agnostic files sit next to this guide. Copy both,
set `COREHPC_PROJECT`, and you should not need to edit either one:

| File | Purpose |
|---|---|
| [`corehpc_nextflow_template.config`](corehpc_nextflow_template.config) | SLURM executor, resource caps, Singularity, reports |
| [`launch_corehpc_template.sh`](launch_corehpc_template.sh) | Submits Nextflow as a driver job; has a built-in smoke test |

Before pointing it at real data, prove the setup works (see section 9):

```bash
export COREHPC_PROJECT=/mnt/gladstone/bioinformatics/projects/<your-project>
./launch_corehpc_template.sh smoke
```

A worked, project-specific version of the same pattern lives in the private
`gladstone-institutes/GB-TH-1359` repo under
`Long_COVID_HIV_Project/step6_process_pmbc_samples_corehpc/`, alongside a
`fetch_results.sh` for pulling the small outputs back off the cluster.

> That launcher pins `NXF_VER` for a parser incompatibility specific to one
> pipeline release. New setups should not copy that line.

---

## 1. How CoreHPC differs from Wynton

These four differences cause essentially every failure below.

| | Wynton | CoreHPC |
|---|---|---|
| Scheduler | SGE (`qsub`) | SLURM (`sbatch`), account `hpc_core`, partition `cpu` |
| Login node | long-running processes tolerated | **kills your processes when SSH drops** (`KillUserProcesses=yes`, no linger) |
| Compute nodes | have internet | **no internet at all** |
| Storage path | `/gladstone/...` | `/mnt/gladstone/...` (same storage, different mount point) |

The path difference matters more than it looks: the Gladstone **data transfer
node** (`gsdt02.gladstone.internal`) mounts the same project as
`/gladstone/...` with no `/mnt` prefix. Config files and rsync commands need
different roots for the same directory.

---

## 2. The driver-job pattern

Because the login node reaps processes on disconnect, you cannot just run
`nextflow run ...` in your SSH session. **Submit Nextflow itself as a SLURM
job.** That "driver" job sits on a compute node orchestrating, and submits one
child SLURM job per pipeline task via the `slurm` executor.

The driver is cheap: 1 CPU, 8 GB, it only polls. Give it a long walltime
(`sinfo -o "%P %l"` shows the partition cap).

```bash
sbatch --parsable \
  --job-name="th1359_rnaseq" \
  --account=hpc_core --partition=cpu \
  --time=5-00:00:00 --mem=8G --cpus-per-task=1 \
  --requeue --signal=B:TERM@300 --open-mode=append \
  --chdir="$SCRIPT_DIR" --output="$LOG" \
  "$JOBSCRIPT"
```

Write the job body to an inspectable `.sbatch` file rather than piping to
`sbatch` on stdin. `launch_corehpc_rnaseq.sh check` writes and prints it
without submitting, which is the fastest way to catch a bad path or an empty
argument before burning a queue slot.

Also print the Nextflow argv from inside the job:

```bash
echo "nextflow argv:"; printf '  <%s>\n' "${NF[@]}"
```

Building the command as a bash array (not backslash continuations) plus this
one line turns "mysterious Nextflow arg error" into "there is an empty string
in position 7".

### Alternative

If the cluster ever permits `loginctl enable-linger $USER`, running the
orchestrator in `tmux` on the login node is simpler, and then only the
short-lived child jobs have walltimes. It was not permitted for us.

---

## 3. Offline execution (the big one)

The driver runs on a compute node with no internet. Three separate things try
to phone home, and each needs its own fix.

**a) The pipeline itself.** `nextflow run nf-core/rnaseq -r 3.19.0` treats that
as a remote project and tries to pull/update it:

```
curl: (7) Failed to connect to www.nextflow.io ...
Pulling nf-core/rnaseq ... Cannot read project manifest -- Cause: No route to host
```

Run the **local** copy produced by `nf-core pipelines download`, with no `-r`:

```bash
nextflow run /path/to/assets/nf-core-rnaseq-3.19.0/3_19_0 ...
```

**b) Self-update and remote configs.** Set `export NXF_OFFLINE=true` in the job
body. This disables the version check and the remote `nf-core/configs` fetch
(nf-core guards that include on `NXF_OFFLINE`), and forces cached
plugins/images only.

**c) Plugins.** `NXF_OFFLINE=true` also means Nextflow will not download the
plugins the pipeline pins, so it fails with `Plugin with id nf-schema not found
in any repository`. Pre-fetch on the login node into the shared home:

```bash
nextflow plugin install nf-schema@2.3.0
ls ~/.nextflow/plugins | grep nf-schema     # expect nf-schema-2.3.0
```

Check the pipeline's `nextflow.config` `plugins { }` block for the exact
versions it pins.

---

## 4. Pre-download containers (login node only)

```bash
export NXF_SINGULARITY_CACHEDIR=$PROJECT/results/singularity-cache

uvx nf-core pipelines download rnaseq \
  --revision 3.19.0 \
  --container-system singularity \
  --container-cache-utilisation amend \
  --compress none \
  --outdir $PROJECT/assets/nf-core-rnaseq-3.19.0
```

The non-obvious flag is **`--container-cache-utilisation amend`**: it puts
images into `$NXF_SINGULARITY_CACHEDIR`, which must equal
`singularity.cacheDir` in your config so the offline run finds them. Confirm
with `ls $NXF_SINGULARITY_CACHEDIR` (expect `.img`/`.sif` files).

Also `nextflow pull nf-core/rnaseq -r 3.19.0` first, to cache the code under
`~/.nextflow` (shared with compute nodes).

---

## 5. The config file

Copy [`corehpc_nextflow_template.config`](corehpc_nextflow_template.config)
rather than writing this from scratch. It reads every path from the environment,
so the only thing you must set is:

```bash
export COREHPC_PROJECT=/mnt/gladstone/bioinformatics/projects/<your-project>
```

Optional overrides: `COREHPC_ACCOUNT`, `COREHPC_PARTITION`, `COREHPC_SCRATCH`,
`COREHPC_WORKDIR`, `COREHPC_CACHEDIR`, `COREHPC_REPORTS`. It throws a clear
error if `COREHPC_PROJECT` is unset, instead of silently writing to the wrong
place.

The parts that matter on CoreHPC, and why:

```groovy
process {
    executor       = 'slurm'
    queue          = partition                    // COREHPC_PARTITION, default 'cpu'
    clusterOptions = "--account=${account}"       // COREHPC_ACCOUNT, default 'hpc_core'

    // Replaces the removed --max_cpus/--max_memory/--max_time params.
    resourceLimits = [ cpus: 4, memory: 60.GB, time: 50.h ]

    // Drop -C (noclobber) from nf-core's default shell. See below.
    shell = ['bash', '-e', '-u', '-o', 'pipefail']
}

singularity {
    enabled    = true
    autoMounts = true
    cacheDir   = cachePath
    runOptions = "--bind ${scratch} --bind ${projRoot}"
}

executor {
    queueSize       = 20        // cap concurrent child jobs
    submitRateLimit = '10 sec'
    exitReadTimeout = '15 min'  // tolerate slow shared-FS exit-file visibility
}

report   { enabled = true; overwrite = true; file = 'reports/execution_report.html' }
timeline { enabled = true; overwrite = true; file = 'reports/execution_timeline.html' }
trace    { enabled = true; overwrite = true; file = 'reports/execution_trace.txt' }
```

### `resourceLimits`, not `--max_cpus`

Recent nf-core pipelines removed `--max_cpus` / `--max_memory` / `--max_time`.
Any request above `resourceLimits` is clamped down instead. Old Wynton launch
commands carrying `--max_*` flags will error.

### noclobber breaks `-resume`

nf-core's default process shell includes `-C` (noclobber). On `-resume` after
an interrupted or requeued run, a task reruns in its **old** work directory
where its output file still exists, and `-C` makes the `>` redirect fail with
`cannot overwrite existing file`. Drop `-C`, keep `-euo pipefail`.

### `overwrite = true` on reports

Without it, a relaunch aborts immediately with `Trace file already exists ...
enable the 'trace.overwrite' option`. Note each relaunch then replaces the
previous run's report files, so copy them aside if you want to keep them.

### Bind mounts

Set both the config `runOptions --bind` and `APPTAINER_BINDPATH` in the job
body, covering scratch and the project root. Containers cannot see paths that
are not bound, and the resulting error usually reads as a missing input file.

### Per-process overrides

Start from the defaults and override only what actually failed. Ours after
tuning:

```groovy
withName: 'STAR_GENOMEGENERATE' { cpus = 4; memory = '60.GB'; time = '8.h' }
withName: 'STAR_ALIGN'          { cpus = 4; memory = '40.GB'; time = '8.h' }

// Loads the decoy-aware index; 16 GB OOM'd (exit 137). Scale on retry.
withName: 'SALMON_QUANT'        { cpus = 4; memory = { 24.GB * task.attempt }; time = '4.h' }

// Fails on upstream fastq md5 mismatches that do not affect alignment.
withName: 'FQ_LINT'             { errorStrategy = 'ignore' }
```

Exit status 137 is always OOM. Prefer `{ N.GB * task.attempt }` over a single
large fixed request, so the first attempt stays schedulable.

---

## 6. Surviving the driver walltime

The driver is a job, so it has a walltime, and a big run can outlast it. The
launcher checkpoints and requeues itself:

- `--signal=B:TERM@300` makes SLURM send `SIGTERM` to the driver 300 s before
  the walltime.
- A trap stops Nextflow gracefully. Nextflow flushes the `-resume` cache and
  cancels its running child jobs.
- `scontrol requeue $SLURM_JOB_ID` restarts the same job, `--requeue` and
  `--open-mode=append` let it resume into the same log, and `-resume` picks up
  from the cache.

```bash
NF_PID=""
on_term() {
  echo "[driver] TERM near walltime -- graceful stop + requeue $(date)"
  kill -TERM "$NF_PID" 2>/dev/null || true
  wait "$NF_PID" 2>/dev/null
  scontrol requeue "$SLURM_JOB_ID" || echo "[driver] requeue failed"
  exit 0
}
trap on_term TERM

"${NF[@]}" &          # background so the trap can fire
NF_PID=$!
wait "$NF_PID"
```

**What survives:** every task that completed and was cached.
**What does not:** tasks mid-flight at the cutoff get recomputed. Nextflow
cannot pause a running task or re-attach to an orphaned child job.

So set `DRIVER_TIME` as high as the partition allows to minimize restarts.
`SIGNAL_GRACE` only needs to cover Nextflow's own shutdown (seconds), not the
running tasks. If the cluster disallows requeue or in-job `scontrol`, the trap
logs a failure and the job ends. Just relaunch manually, `-resume` still works.

Always run with `-resume`. It costs nothing on a fresh run and saves days on a
restart.

---

## 7. Samplesheet gotcha: sample IDs must not be bare numbers

Sample IDs like `4107` are parsed by nf-schema as integers and fail the
samplesheet's `type: string` check. Prefix them (`TH_4107`). This also keeps R
from mangling counts-matrix column names later.

---

## 8. Getting results back

Transfer through **`gsdt02.gladstone.internal`**, not a CoreHPC login node, and
remember the mount point changes to `/gladstone/...` there.

Pull only the small outputs. Alignments, fastqs, bigwigs and the on-the-fly
STAR/RSEM/salmon index are tens to hundreds of GB and belong on the cluster.
`fetch_results.sh` excludes them by pattern plus a `--max-size` backstop, and
its `check` mode prints exactly what the size cap would skip so nothing
disappears silently.

Note that the Nextflow `report`/`timeline`/`trace` files land next to the
launcher (per the `file =` paths in the config), **not** in `--outdir`, so
they need a second rsync.

---

## 9. Verifying the setup

Do not debug your cluster setup and your pipeline at the same time. Run the
smoke test first:

```bash
export COREHPC_PROJECT=/mnt/gladstone/bioinformatics/projects/<your-project>
./launch_corehpc_template.sh smoke check   # inspect the sbatch, submit nothing
./launch_corehpc_template.sh smoke         # submit it
```

It writes a three-task pipeline into `$COREHPC_PROJECT/corehpc_smoke_test/` and
runs it through the same driver job, config and SLURM executor a real run uses.
It takes a couple of minutes. Then:

```bash
cat $COREHPC_PROJECT/corehpc_smoke_test/results/smoke_report.txt
```

Every task should report a real SLURM job id, the cpus and memory your config
granted it, and a work directory under `$COREHPC_PROJECT`. That confirms the
scheduler, account, partition, config and work directory are all wired up. If
the container cache already holds an image, the report also confirms Singularity
runs and that the project bind mount is visible from inside the container.

Now run it a **second** time. Every task should say `CACHED`, which is what
proves `-resume` works. Without that, a requeue after the driver walltime will
recompute everything.

Delete `$COREHPC_PROJECT/corehpc_smoke_test/` when you are done.

### Pre-flight checklist

Once the smoke test passes, run this before each new pipeline:

```bash
which nextflow; nextflow -version                   # engine present
which singularity || which apptainer                # if only apptainer, use -profile apptainer
ls ~/.nextflow/plugins                              # pinned plugins pre-fetched
ls $NXF_SINGULARITY_CACHEDIR                        # images present (.img/.sif)
ls $PIPELINE/main.nf                                # local pipeline copy exists
ls $SAMPLESHEET                                     # and every path inside it
sinfo -o "%P %l"                                    # partition walltime cap vs DRIVER_TIME
sacctmgr show assoc user=$USER format=account,partition   # account is valid
```

Then run `./launch_corehpc_template.sh check`, read the printed sbatch, and submit.

Watch it with:

```bash
squeue -u $USER            # driver first, then child nf-* jobs
tail -f <driver log>
```

---

## 10. Error quick reference

| Symptom | Cause | Fix |
|---|---|---|
| `Cannot read project manifest -- No route to host` | Compute node pulling the pipeline from GitHub | Run the local downloaded copy, no `-r`; `NXF_OFFLINE=true` |
| `Plugin with id nf-schema not found in any repository` | Offline mode cannot download plugins | `nextflow plugin install nf-schema@2.3.0` on the login node |
| Run dies when SSH disconnects | Login node reaps user processes | Submit Nextflow as a driver job |
| `cannot overwrite existing file` on `-resume` | nf-core default shell includes `-C` (noclobber) | `shell = ['bash','-e','-u','-o','pipefail']` |
| `Trace file already exists` | Report files from the previous run | `overwrite = true` in `report`/`timeline`/`trace` |
| Exit status 137 | OOM | Raise memory, ideally `{ N.GB * task.attempt }` |
| `--max_cpus` unrecognized | Params removed from recent nf-core | `resourceLimits = [...]` in the config |
| Samplesheet schema error on a sample ID | Bare-numeric ID read as an integer | Prefix the ID (`TH_4107`) |
| Input file "not found" inside a container | Path not bound | Add to `runOptions --bind` and `APPTAINER_BINDPATH` |
| Job dies at driver walltime | Run outlasted `DRIVER_TIME` | Requeue trap above, and relaunch with `-resume` |

---

## 11. Porting a Wynton launch command to CoreHPC

1. `/gladstone/...` becomes `/mnt/gladstone/...` (but `/gladstone/...` again on the transfer node).
2. `qsub` config becomes `executor = 'slurm'` + `queue` + `clusterOptions = '--account=hpc_core'`.
3. Drop `--max_cpus/--max_memory/--max_time`, add `resourceLimits`.
4. Drop `module load openjdk/17`; Java is already on PATH.
5. Wrap the whole thing in a driver job instead of running it in the shell.
6. Pre-download pipeline, plugins and containers on the login node.
7. Add `NXF_OFFLINE=true` and point at the local pipeline copy.
