#!/bin/bash
# Runs all the numerical experiments sequentially, one MATLAB
# process per experiment and size (so that, if one run fails or runs out of
# memory, the others are not affected).
#
# Output, in results/<machine>_<date>/ (or $OUT_DIR if set):
#   summary.md     comparison of compress, two-pass and restart, results of all the methods, exit status of every run
#   <test>.csv     results of every experiment (all sizes)
#   raw/           results saved by the test scripts after every method (one file per size)
#   logs/          MATLAB output of every run
#   run.log        progress of the runs
#   jobs.csv       exit status and duration of every run
#   machine.txt    machine, MATLAB version and settings
#
# Environment variables:
#   MATLAB     MATLAB executable (default: matlab)
#   NTHREADS   number of computational threads used by MATLAB (default: 4)
#   MEMLIMIT   memory limit for every MATLAB run, e.g. 12G (default: none);
#              requires systemd-run, a run exceeding it is killed (exit status 137)
#   OUT_DIR    output folder (default: results/<machine>_<date>)
#   JOBS       experiments to run, as script:number_of_sizes (default: all)

cd "$(dirname "$0")" || exit 1
ROOT=$(pwd)
MATLAB=${MATLAB:-matlab}
NTHREADS=${NTHREADS:-4}
OUT_DIR=${OUT_DIR:-results/$(hostname -s)_$(date +%F)}
mkdir -p "$OUT_DIR/raw" "$OUT_DIR/logs"
OUT_DIR=$(cd "$OUT_DIR" && pwd)
export LYAP_RESULTS_DIR="$OUT_DIR/raw"
exec > >(tee -a "$OUT_DIR/run.log") 2>&1   # progress of the runs, also saved in run.log

PREFIX=()
if [ -n "$MEMLIMIT" ]; then
    PREFIX=(systemd-run --user --scope --quiet -p MemoryMax="$MEMLIMIT" -p MemorySwapMax=0)
fi

{
    echo "date:      $(date '+%F %T')"
    echo "host:      $(hostname -f)"
    echo "cpu:       $(lscpu | grep 'Model name' | sed 's/.*: *//')"
    echo "cores:     $(nproc)"
    echo "memory:    $(free -g | awk '/Mem:/ {print $2 " GB"}')"
    echo "threads:   $NTHREADS (maxNumCompThreads)"
    echo "mem limit: ${MEMLIMIT:-none}"
    echo "matlab:    $("$MATLAB" -batch "disp(version)" 2>/dev/null | tail -1)"
} > "$OUT_DIR/machine.txt"

echo "script,size,status,seconds" > "$OUT_DIR/jobs.csv"

# script:number of sizes (0 = script without sizes)
JOBS=${JOBS:-"test_2DLap:4 test_2DLap_sensitivity:0 test_3DLap:4 test_fenics_rail:3 test_heat:3"}

for job in $JOBS; do
    t=${job%%:*}
    nsizes=${job##*:}
    if [ "$nsizes" -eq 0 ]; then sizes="all"; else sizes=$(seq 1 "$nsizes"); fi
    for i in $sizes; do
        if [ "$i" = "all" ]; then cmd="$t"; log="$t"; else cmd="size_idx = $i; $t"; log="${t}_size$i"; fi
        echo "[$(date '+%F %T')] start $log"
        start=$(date +%s)
        "${PREFIX[@]}" "$MATLAB" -batch "maxNumCompThreads($NTHREADS); $cmd" > "$OUT_DIR/logs/$log.log" 2>&1
        status=$?
        echo "$t,$i,$status,$(( $(date +%s) - start ))" >> "$OUT_DIR/jobs.csv"
        echo "[$(date '+%F %T')] end   $log (exit status $status)"
    done
done

"$MATLAB" -batch "addpath('$ROOT'); setup_paths; make_summary('$OUT_DIR')" > "$OUT_DIR/logs/make_summary.log" 2>&1
echo "[$(date '+%F %T')] all experiments finished, see $OUT_DIR/summary.md"
