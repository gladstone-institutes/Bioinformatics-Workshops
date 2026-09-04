#!/bin/bash
# ============================================================================
# Nextflow driver-job launcher for UCSF CoreHPC (SLURM). Project-agnostic.
#
# The CoreHPC login node kills your processes when the SSH session drops, so
# Nextflow itself is submitted as a batch job ("the driver"). It sits on a
# compute node and submits one child job per task via the slurm executor.
#
# Usage:
#   export COREHPC_PROJECT=/mnt/gladstone/bioinformatics/projects/<your-project>
#
#   ./launch_corehpc_template.sh smoke   # verify the setup end to end, no real data
#   ./launch_corehpc_template.sh check   # write + print the driver sbatch, submit nothing
#   ./launch_corehpc_template.sh smoke check   # same, for the smoke test
#   ./launch_corehpc_template.sh         # submit the real run
#
# Real runs need PIPELINE and SAMPLESHEET, and usually EXTRA:
#   PIPELINE=$COREHPC_PROJECT/assets/nf-core-rnaseq-3.19.0/3_19_0 \
#   SAMPLESHEET=$COREHPC_PROJECT/input/samplesheet.csv \
#   EXTRA="--fasta $REF/genome.fa --gtf $REF/genes.gtf.gz" \
#     ./launch_corehpc_template.sh
#
# See running_nextflow_on_corehpc.md for the login-node prep this assumes.
# ============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

MODE="${1:-run}"

: "${COREHPC_PROJECT:?Set COREHPC_PROJECT to your project root (note the /mnt prefix)}"
[ -d "$COREHPC_PROJECT" ] || { echo "Error: COREHPC_PROJECT does not exist: $COREHPC_PROJECT" >&2; exit 1; }

CONFIG="${CONFIG:-$SCRIPT_DIR/corehpc_nextflow_template.config}"
[ -f "$CONFIG" ] || { echo "Error: config not found: $CONFIG" >&2; exit 1; }

SCRATCH="${COREHPC_SCRATCH:-/mnt/scratch}"
CACHEDIR="${COREHPC_CACHEDIR:-$COREHPC_PROJECT/singularity-cache}"
WORKDIR="${COREHPC_WORKDIR:-$COREHPC_PROJECT/work}"
PROFILE="${PROFILE:-singularity}"

# Driver job knobs.
ACCOUNT="${COREHPC_ACCOUNT:-hpc_core}"
PARTITION="${COREHPC_PARTITION:-cpu}"
DRIVER_TIME="${DRIVER_TIME:-5-00:00:00}"   # must outlast the whole run
DRIVER_MEM="${DRIVER_MEM:-8G}"             # the driver only polls
SIGNAL_GRACE="${SIGNAL_GRACE:-300}"        # secs before walltime to stop + requeue

# JVM heap for the driver, derived from DRIVER_MEM so the two cannot drift apart.
# The JVM sizes its default heap from the node's RAM, not the SLURM cgroup, so an
# uncapped heap grows past --mem and SLURM kills the driver with exit 137. Leave
# ~25% for non-heap overhead (metaspace, thread stacks, GC, code cache).
_mem_gb="${DRIVER_MEM%[Bb]}"; _mem_gb="${_mem_gb%[Gg]}"
if [[ "$_mem_gb" =~ ^[0-9]+$ ]] && (( _mem_gb >= 2 )); then
  NXF_HEAP="${NXF_HEAP:-$(( _mem_gb * 3 / 4 ))g}"
else
  NXF_HEAP="${NXF_HEAP:-4g}"
  echo "Note: could not parse DRIVER_MEM=$DRIVER_MEM, defaulting the JVM heap to $NXF_HEAP." >&2
fi

require() { command -v "$1" &>/dev/null || { echo "Error: $1 not found." >&2; exit 1; }; }
require nextflow

# --- smoke mode: a tiny pipeline that exercises the same path as a real run ---
if [[ "$MODE" == "smoke" ]]; then
  SMOKE_DIR="$COREHPC_PROJECT/corehpc_smoke_test"
  mkdir -p "$SMOKE_DIR"
  PIPELINE="$SMOKE_DIR/smoke.nf"
  OUTDIR="${OUTDIR:-$SMOKE_DIR/results}"
  WORKDIR="$SMOKE_DIR/work"
  DRIVER_TIME="00:30:00"

  # Use any image already in the cache to prove singularity + binds work. If the
  # cache is empty the container check is skipped, and the rest still runs.
  SMOKE_IMG="$(ls "$CACHEDIR"/*.img "$CACHEDIR"/*.sif 2>/dev/null | head -1 || true)"

  cat > "$PIPELINE" <<'NF'
// Smoke test for a CoreHPC Nextflow setup. Submits real SLURM jobs, writes a
// report, and caches so a second run can be checked with -resume.
nextflow.enable.dsl = 2

params.outdir    = 'results'
params.container = ''
params.bindcheck = ''

process PING {
    tag "task_${idx}"
    input:
    val idx
    output:
    path "ping_${idx}.txt"
    script:
    """
    { echo "task        : ${idx}"
      echo "host        : \$(hostname)"
      echo "slurm job   : \${SLURM_JOB_ID:-NOT-A-SLURM-JOB}"
      echo "cpus given  : ${task.cpus}"
      echo "memory given: ${task.memory}"
      echo "workdir     : \$(pwd)"
    } > ping_${idx}.txt
    """
}

process CONTAINER_CHECK {
    tag "container"
    container params.container
    output:
    path 'container.txt'
    when:
    params.container != ''
    script:
    """
    if ls ${params.bindcheck} > /dev/null 2>&1; then
        BIND="OK, can see ${params.bindcheck}"
    else
        BIND="FAILED, cannot see ${params.bindcheck} from inside the container"
    fi
    { echo "image       : ${params.container}"
      echo "container id: \${SINGULARITY_NAME:-\${APPTAINER_NAME:-none, ran on the host}}"
      echo "bind mount  : \$BIND"
    } > container.txt
    """
}

process REPORT {
    tag "report"
    publishDir params.outdir, mode: 'copy', overwrite: true
    input:
    path parts
    output:
    path 'smoke_report.txt'
    script:
    """
    { echo "CoreHPC Nextflow smoke test"
      echo "generated: \$(date)"
      echo "==========================================="
      cat ${parts}
      echo "==========================================="
      echo "If every task above shows a real SLURM job id, the scheduler,"
      echo "config and work directory are wired up correctly."
    } > smoke_report.txt
    """
}

workflow {
    pings = PING( Channel.of(1, 2, 3) )
    REPORT( pings.mix( CONTAINER_CHECK() ).collect() )
}
NF

  EXTRA=""
  if [[ -n "$SMOKE_IMG" ]]; then
    EXTRA="--container $SMOKE_IMG --bindcheck $COREHPC_PROJECT"
    echo "Smoke test will check containers with: $SMOKE_IMG"
  else
    echo "Note: no .img/.sif in $CACHEDIR, so the container check is skipped."
    echo "      Scheduler, config and work dir are still exercised."
  fi
  SAMPLESHEET=""
  echo "Smoke test dir: $SMOKE_DIR"
  echo ""
fi

# --- real runs ---------------------------------------------------------------
PIPELINE="${PIPELINE:-}"
SAMPLESHEET="${SAMPLESHEET-}"
OUTDIR="${OUTDIR:-$COREHPC_PROJECT/results}"
EXTRA="${EXTRA:-}"

[ -n "$PIPELINE" ] || { echo "Error: set PIPELINE to your local pipeline copy (see README)." >&2; exit 1; }
[ -f "$PIPELINE" ] || [ -f "$PIPELINE/main.nf" ] || {
  echo "Error: no pipeline at $PIPELINE (expected a .nf file or a dir with main.nf)." >&2
  echo "  Compute nodes have no internet, so this must be a LOCAL copy from" >&2
  echo "  'nf-core pipelines download', not a remote name like nf-core/rnaseq." >&2
  exit 1; }

if [[ -n "$SAMPLESHEET" ]]; then
  [ -f "$SAMPLESHEET" ] || { echo "Error: samplesheet not found: $SAMPLESHEET" >&2; exit 1; }
  EXTRA="--input $SAMPLESHEET $EXTRA"
fi

DRIVERDIR="$OUTDIR/logs/driver"
mkdir -p "$OUTDIR" "$WORKDIR" "$CACHEDIR" "$DRIVERDIR"

STAMP=$(date +%Y%m%d_%H%M%S)
JOBSCRIPT="$DRIVERDIR/nf_${STAMP}.sbatch"
LOG="$DRIVERDIR/nf_${STAMP}_%j.log"

# Write the driver body to a file so it can be inspected before submission.
cat > "$JOBSCRIPT" <<EOF
#!/bin/bash
set -uo pipefail
cd "$SCRIPT_DIR"
echo "[driver] \$SLURM_JOB_ID on \$(hostname), started \$(date)"

export NXF_OPTS='-Xms1g -Xmx$NXF_HEAP'
export COREHPC_PROJECT="$COREHPC_PROJECT"
export COREHPC_SCRATCH="$SCRATCH"
export COREHPC_ACCOUNT="$ACCOUNT"
export COREHPC_PARTITION="$PARTITION"
export COREHPC_WORKDIR="$WORKDIR"
export COREHPC_CACHEDIR="$CACHEDIR"
export NXF_SINGULARITY_CACHEDIR="$CACHEDIR"
export APPTAINER_BINDPATH="$SCRATCH,$COREHPC_PROJECT"
# Compute nodes have no internet: no self-update check, no remote nf-core/configs
# fetch, cached plugins and images only.
export NXF_OFFLINE=true

# Walltime survival. SLURM sends SIGTERM $SIGNAL_GRACE secs before the walltime;
# stop Nextflow gracefully so it flushes the -resume cache and cancels its child
# jobs, then requeue this job to continue. Completed tasks are kept; tasks that
# were mid-flight are recomputed.
NF_PID=""
on_term() {
  echo "[driver] TERM near walltime, graceful stop + requeue \$(date)"
  kill -TERM "\$NF_PID" 2>/dev/null || true
  wait "\$NF_PID" 2>/dev/null
  scontrol requeue "\$SLURM_JOB_ID" || echo "[driver] requeue failed (is --requeue allowed here?)"
  exit 0
}
trap on_term TERM

NF=(nextflow run "$PIPELINE" -profile $PROFILE \
  --outdir "$OUTDIR" \
  -work-dir "$WORKDIR" \
  -c "$CONFIG" \
  -resume \
  $EXTRA)

# Print argv so a stray or empty argument is visible in the log.
echo "nextflow argv:"
printf '  <%s>\n' "\${NF[@]}"

# Background it so the trap can act on the walltime signal.
"\${NF[@]}" &
NF_PID=\$!
wait "\$NF_PID"
rc=\$?
echo "[driver] Nextflow exited rc=\$rc \$(date)"
exit \$rc
EOF

if [[ "$MODE" == "check" || "${2:-}" == "check" ]]; then
  echo "Driver job script written (NOT submitted): $JOBSCRIPT"
  echo "----------------------------------------------------------------------"
  cat "$JOBSCRIPT"
  exit 0
fi

require sbatch

JOBID=$(sbatch --parsable \
  --job-name="${JOB_NAME:-nf_driver}" \
  --account="$ACCOUNT" \
  --partition="$PARTITION" \
  --time="$DRIVER_TIME" \
  --mem="$DRIVER_MEM" \
  --cpus-per-task=1 \
  --requeue \
  --signal="B:TERM@$SIGNAL_GRACE" \
  --open-mode=append \
  --chdir="$SCRIPT_DIR" \
  --output="$LOG" \
  "$JOBSCRIPT")

RESOLVED_LOG="$DRIVERDIR/nf_${STAMP}_${JOBID}.log"
echo "Submitted driver job $JOBID"
echo "  driver:   $DRIVER_TIME walltime, $DRIVER_MEM (JVM heap $NXF_HEAP), 1 cpu, partition $PARTITION, account $ACCOUNT"
echo "  pipeline: $PIPELINE"
echo "  outdir:   $OUTDIR"
echo "  log:      $RESOLVED_LOG"
echo ""
echo "Watch:"
echo "  squeue -u \$USER          # driver first, then child nf-* jobs"
echo "  tail -f $RESOLVED_LOG"

if [[ "$MODE" == "smoke" ]]; then
  echo ""
  echo "When it finishes, read:"
  echo "  cat $OUTDIR/smoke_report.txt"
  echo "Then run it again to confirm caching works. Every task should say CACHED."
fi
