#!/usr/bin/env bash
#
# Finds the largest batch size that fits in GPU memory by running your script
# with different sizes (binary search) until it stops OOMing. Run inside a GPU
# salloc session.
#
# IMPORTANT: your script must read the size from the BATCH_SIZE env var and 
# exit 0 on success:
#   import os
#   batch_size = int(os.environ.get("BATCH_SIZE", 32))
#
# Usage: edit the Settings block below, then run: ./find_batch_size.sh

# ----------------------------------------------------------------------------
# Settings: change these values, then save and run the script.
# ----------------------------------------------------------------------------
PYTHON_SCRIPT=train.py   # name of your training script
START_BATCH_SIZE=1024    # largest size to try; the result is capped here, so set it above any size that might fit
MIN_BATCH_SIZE=8         # smallest size to try; if even this runs out of memory, nothing fit
STEP=8                   # only test multiples of this number
# ----------------------------------------------------------------------------

if [ ! -f "$PYTHON_SCRIPT" ]; then
  echo "Error: Python script '$PYTHON_SCRIPT' not found."
  echo "Set PYTHON_SCRIPT in the Settings block at the top of this file."
  exit 1
fi

# Snap bounds to multiples of STEP.
START_BATCH_SIZE=$(( START_BATCH_SIZE / STEP * STEP ))
MIN_BATCH_SIZE=$(( MIN_BATCH_SIZE / STEP * STEP ))
if [ "$MIN_BATCH_SIZE" -lt "$STEP" ]; then MIN_BATCH_SIZE=$STEP; fi

OOM_PATTERN="out of memory|CUDA error: out of memory|OutOfMemoryError|CUBLAS_STATUS_ALLOC_FAILED"
RUN_LOG=$(mktemp)
trap 'rm -f "$RUN_LOG"' EXIT

# Binary search: assumes if a size fits, every smaller size fits too.
low=$MIN_BATCH_SIZE
high=$START_BATCH_SIZE
best=-1

echo "Searching $MIN_BATCH_SIZE to $START_BATCH_SIZE (multiples of $STEP), script: $PYTHON_SCRIPT"
echo

while [ "$low" -le "$high" ]; do
  candidate=$(( (low + high) / 2 ))
  candidate=$(( candidate / STEP * STEP ))
  if [ "$candidate" -lt "$MIN_BATCH_SIZE" ]; then candidate=$MIN_BATCH_SIZE; fi

  echo ">> Trying batch size $candidate"

  # tee shows output live and saves it for the OOM check; PIPESTATUS[0] is python's code.
  BATCH_SIZE="$candidate" python "$PYTHON_SCRIPT" 2>&1 | tee "$RUN_LOG"
  exit_code=${PIPESTATUS[0]}

  if [ "$exit_code" -eq 0 ]; then
    echo ">> OK: $candidate fits, trying larger."
    best=$candidate
    low=$(( candidate + STEP ))
  elif [ "$exit_code" -eq 137 ] || grep -qiE "$OOM_PATTERN" "$RUN_LOG"; then
    echo ">> OOM at $candidate, trying smaller."
    high=$(( candidate - STEP ))
  else
    echo "Error: failed at $candidate for a non-OOM reason (exit $exit_code). Fix it and re-run."
    exit "$exit_code"
  fi
  echo
done

if [ "$best" -ge 0 ]; then
  echo "Largest batch size that fits: $best"
else
  echo "Nothing fit: even $MIN_BATCH_SIZE ran out of memory."
fi
