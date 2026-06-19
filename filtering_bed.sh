#!/usr/bin/env bash

# ==============================================================================
# Script Name: filter_genomic_regions.sh
# Description: Automated pipeline filtering genomic target regions (WES or 
#              targeted panels). The script detects the input format (BED or 
#              Picard interval_list), standardizes coordinates to 0-based BED, 
#              and subtracts regions associated with Simple Motifs and 
#              Homopolymers to minimize sequencing artifacts.
# Dependencies: bedtools (v2.25.0+), awk
# ==============================================================================

set -euo pipefail

# --- Configuration ---
# Input file (can be .bed or .interval_list)
INPUT_TARGETS="path/to/your_file.interval_list"

# Homopolymers and Simple Repeats BED files
BED_SM="path/to/rmsk_SimpleMotifs_only_slop1.bed"
BED_HOMO="path/to/GRCh38_AllHomopolymers_gt6bp_imperfectgt10bp_slop1.bed"

# Output
OUT_DIR="./filtered_exome_intervals"
mkdir -p "$OUT_DIR"

EXOME_BED="$OUT_DIR/input_standardized.bed"
INTERMEDIATE_BED="$OUT_DIR/targets_noSimpleMotifs.bed"
FINAL_BED="$OUT_DIR/targets_filtered_final.bed"

echo "-------------------------------------------------------"
echo "Starting Genomic Region Filtering Pipeline"
echo "Timestamp: $(date)"
echo "-------------------------------------------------------"

# --- Step 1: Format Detection and Normalization ---
echo "[Step 1/3] Checking input format..."

if [[ "$INPUT_TARGETS" == *.interval_list ]]; then
    echo "Detected Picard interval_list. Converting to 0-based BED..."
    grep -v '^@' "$INPUT_TARGETS" | awk -v OFS='\t' '{print $1, $2-1, $3}' > "$EXOME_BED"
elif [[ "$INPUT_TARGETS" == *.bed ]]; then
    echo "Detected BED file. Proceeding with copy..."
    cp "$INPUT_TARGETS" "$EXOME_BED"
else
    echo "Error: Input file must end in .bed or .interval_list"
    exit 1
fi

# --- Step 2: Remove Simple Motifs ---
echo "[Step 2/3] Subtracting Simple Motifs..."
bedtools subtract -a "$EXOME_BED" -b "$BED_SM" > "$INTERMEDIATE_BED"

# --- Step 3: Remove Homopolymers ---
echo "[Step 3/3] Subtracting Homopolymers..."
bedtools subtract -a "$INTERMEDIATE_BED" -b "$BED_HOMO" > "$FINAL_BED"

# --- Summary Statistics ---
echo "-------------------------------------------------------"
echo "Filtering Summary (Total Intervals):"
echo "  Initial regions:      $(wc -l < "$EXOME_BED")"
echo "  Final Filtered:       $(wc -l < "$FINAL_BED")"
echo "-------------------------------------------------------"
echo "Process completed. Output: $FINAL_BED"