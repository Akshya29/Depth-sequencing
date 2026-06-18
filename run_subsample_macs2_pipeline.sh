#!/usr/bin/env bash
#
# run_subsample_macs2_pipeline.sh
#
# Subsample sorted BAM files at a range of percentages, call peaks on each
# subsample with MACS2, and tabulate peak counts vs. subsampling percentage.
#
# Usage:
#   ./run_subsample_macs2_pipeline.sh -i /path/to/bam_dir [options]
#
# Required:
#   -i DIR        Directory containing *_sorted.bam input files
#
# Optional:
#   -o DIR        Output directory (default: ./subsample_macs2_output)
#   -s START      Start percentage (default: 5)
#   -e END        End percentage, inclusive (default: 100)
#   -t STEP       Step size in percent (default: 15)
#   -g GENOME     MACS2 genome size, e.g. 2.3e9, hs, mm (default: 2.3e9)
#   -d SEED       Subsampling seed (default: 42)
#   -m MODULE     Environment module to load for samtools, e.g.
#                 samtools/1.18-gcc-12.3.0 (optional; skipped if omitted
#                 or if `module` is not available)
#   -j JOBS       Number of BAMs to process in parallel (default: 1)
#   -h            Show this help message and exit
#
# Requires: samtools, macs2, awk, bc (or awk-only fallback), GNU grep (-P)
#
# Example:
#   ./run_subsample_macs2_pipeline.sh \
#       -i /QRISdata/Q3338/Akshya/UMR-seq-aim2-maize-deep/analysis/trimmed_align_bowtie2 \
#       -o /scratch/akshya/subsample_run1 \
#       -s 5 -e 100 -t 15 \
#       -g 2.3e9 \
#       -m samtools/1.18-gcc-12.3.0 \
#       -j 4
#
set -uo pipefail

# ---------------------------- defaults -------------------------------------
IN_DIR=""
OUT_DIR="./subsample_macs2_output"
START_PCT=5
END_PCT=100
STEP_PCT=15
GENOME_SIZE="2.3e9"
SEED=42
MODULE_NAME=""
JOBS=1

# ---------------------------- usage -----------------------------------------
usage() {
    sed -n '2,33p' "$0" | sed 's/^# \{0,1\}//'
    exit "${1:-0}"
}

# ---------------------------- parse args -------------------------------------
while getopts ":i:o:s:e:t:g:d:m:j:h" opt; do
    case "$opt" in
        i) IN_DIR="$OPTARG" ;;
        o) OUT_DIR="$OPTARG" ;;
        s) START_PCT="$OPTARG" ;;
        e) END_PCT="$OPTARG" ;;
        t) STEP_PCT="$OPTARG" ;;
        g) GENOME_SIZE="$OPTARG" ;;
        d) SEED="$OPTARG" ;;
        m) MODULE_NAME="$OPTARG" ;;
        j) JOBS="$OPTARG" ;;
        h) usage 0 ;;
        \?) echo "Unknown option: -$OPTARG" >&2; usage 1 ;;
        :) echo "Option -$OPTARG requires an argument" >&2; usage 1 ;;
    esac
done

# ---------------------------- validation --------------------------------------
if [[ -z "$IN_DIR" ]]; then
    echo "ERROR: input directory (-i) is required." >&2
    usage 1
fi
if [[ ! -d "$IN_DIR" ]]; then
    echo "ERROR: input directory does not exist: $IN_DIR" >&2
    exit 1
fi

for tool in samtools macs2 awk; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        echo "ERROR: required tool '$tool' not found in PATH." >&2
        echo "       Load it via a module (-m) or your environment before running." >&2
        exit 1
    fi
done

if [[ "$START_PCT" -lt 1 || "$END_PCT" -gt 100 || "$START_PCT" -gt "$END_PCT" ]]; then
    echo "ERROR: percentage range invalid (start=$START_PCT, end=$END_PCT). Must be 1-100, start<=end." >&2
    exit 1
fi

# ---------------------------- optional module load -----------------------------
if [[ -n "$MODULE_NAME" ]]; then
    if command -v module >/dev/null 2>&1; then
        echo "Loading module: $MODULE_NAME"
        # shellcheck disable=SC1090
        module load "$MODULE_NAME" || { echo "ERROR: failed to load module $MODULE_NAME" >&2; exit 1; }
    else
        echo "WARNING: 'module' command not available; skipping module load for $MODULE_NAME" >&2
    fi
fi

# ---------------------------- set up output dirs --------------------------------
SUBSAMPLE_DIR="${OUT_DIR}/subsampled_bams"
PEAKS_DIR="${OUT_DIR}/macs2_peaks"
mkdir -p "$SUBSAMPLE_DIR" "$PEAKS_DIR"

LOG_FILE="${OUT_DIR}/pipeline.log"
echo "Pipeline started: $(date)" | tee "$LOG_FILE"
echo "Input dir:        $IN_DIR" | tee -a "$LOG_FILE"
echo "Output dir:       $OUT_DIR" | tee -a "$LOG_FILE"
echo "Percent range:    ${START_PCT}-${END_PCT} step ${STEP_PCT}" | tee -a "$LOG_FILE"
echo "Genome size:      $GENOME_SIZE" | tee -a "$LOG_FILE"
echo "Seed:             $SEED" | tee -a "$LOG_FILE"
echo "Parallel jobs:    $JOBS" | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"

# ---------------------------- gather input BAMs ----------------------------------
shopt -s nullglob
mapfile -t BAM_FILES < <(find "$IN_DIR" -maxdepth 1 -name "*_sorted.bam" ! -name "*subsampled*")

if [[ "${#BAM_FILES[@]}" -eq 0 ]]; then
    echo "ERROR: no *_sorted.bam files found in $IN_DIR (excluding files containing 'subsampled')." >&2
    exit 1
fi

echo "Found ${#BAM_FILES[@]} input BAM(s):" | tee -a "$LOG_FILE"
printf '  %s\n' "${BAM_FILES[@]}" | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"

# ---------------------------- worker function -------------------------------------
process_one_bam() {
    local bam="$1"
    local base
    base=$(basename "$bam" .bam)

    local pct
    for pct in $(seq "$START_PCT" "$STEP_PCT" "$END_PCT"); do
        local fraction
        fraction=$(awk -v p="$pct" 'BEGIN { printf "%.2f", p/100 }')

        local out_file="${SUBSAMPLE_DIR}/${base}_subsampled_${pct}p.bam"

        echo "[$(date '+%H:%M:%S')] Subsampling $base at ${pct}% -> $(basename "$out_file")"

        if ! samtools view -b -s "${SEED}.${fraction#*.}" "$bam" > "$out_file"; then
            echo "ERROR: samtools view failed for $bam at ${pct}%" >&2
            rm -f "$out_file"
            continue
        fi
        if ! samtools index "$out_file"; then
            echo "ERROR: samtools index failed for $out_file" >&2
            continue
        fi

        local name
        name=$(basename "$out_file" .bam)

        echo "[$(date '+%H:%M:%S')] Calling MACS2 peaks for $name"
        macs2 callpeak \
            -t "$out_file" \
            -f BAM \
            -g "$GENOME_SIZE" \
            -n "$name" \
            --outdir "$PEAKS_DIR" \
            >> "${PEAKS_DIR}/${name}.macs2.log" 2>&1

        if [[ $? -ne 0 ]]; then
            echo "ERROR: macs2 callpeak failed for $name (see ${PEAKS_DIR}/${name}.macs2.log)" >&2
        fi
    done
}

export -f process_one_bam
export SUBSAMPLE_DIR PEAKS_DIR SEED GENOME_SIZE START_PCT STEP_PCT END_PCT

# ---------------------------- run (parallel or serial) ------------------------------
if [[ "$JOBS" -gt 1 ]] && command -v parallel >/dev/null 2>&1; then
    echo "Running with GNU parallel ($JOBS jobs)..." | tee -a "$LOG_FILE"
    printf '%s\n' "${BAM_FILES[@]}" | parallel -j "$JOBS" process_one_bam {} 2>&1 | tee -a "$LOG_FILE"
else
    if [[ "$JOBS" -gt 1 ]]; then
        echo "NOTE: GNU 'parallel' not found; falling back to serial processing." | tee -a "$LOG_FILE"
    fi
    for bam in "${BAM_FILES[@]}"; do
        process_one_bam "$bam" 2>&1 | tee -a "$LOG_FILE"
    done
fi

# ---------------------------- tabulate peak counts -----------------------------------
echo "" | tee -a "$LOG_FILE"
echo "Tabulating peak counts..." | tee -a "$LOG_FILE"

PEAK_COUNTS_FILE="${OUT_DIR}/peak_counts.tsv"
{
    echo -e "percent\tsample\tpeak_count"
    for f in "$PEAKS_DIR"/*_peaks.narrowPeak; do
        [[ -e "$f" ]] || continue
        base=$(basename "$f" "_peaks.narrowPeak")
        pct=$(echo "$base" | grep -oP '(?<=_subsampled_)\d+(?=p)')
        count=$(wc -l < "$f")
        echo -e "${pct}\t${base}\t${count}"
    done | sort -n
} > "$PEAK_COUNTS_FILE"

echo "" | tee -a "$LOG_FILE"
echo "Done: $(date)" | tee -a "$LOG_FILE"
echo "Peak counts written to: $PEAK_COUNTS_FILE" | tee -a "$LOG_FILE"
echo "Subsampled BAMs in:     $SUBSAMPLE_DIR" | tee -a "$LOG_FILE"
echo "MACS2 outputs in:       $PEAKS_DIR" | tee -a "$LOG_FILE"
