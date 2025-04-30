#!/bin/bash

# Exit on error
set -e

# Input parameters
CSV_FILE="$1"
DEMUX_DIR="$2"
OUTPUT_DIR="$3"

# Check arguments
if [[ -z "$CSV_FILE" || -z "$DEMUX_DIR" || -z "$OUTPUT_DIR" ]]; then
  echo "Usage: $0 <CSV_FILE> <DEMUX_DIR> <OUTPUT_DIR>"
  exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Temp file to track species/sample mappings
TMP_FILE=$(mktemp)

echo "$CSV_FILE"

# Read the CSV and build the list
while IFS=',' read -r species sample; do
  sample=$(echo "$sample" | tr -d '\r')
  species=$(echo "$species" | tr -d '\r')
  echo "$species,$DEMUX_DIR/${sample}_R1.fastq.gz,$DEMUX_DIR/${sample}_R2.fastq.gz" >> "$TMP_FILE"
done < "$CSV_FILE"

# Get list of unique species
species_list=$(cut -d',' -f1 "$TMP_FILE" | sort | uniq)

# Process each species
for species in $species_list; do
  echo "Processing $species..."

  R1_FILES=()
  R2_FILES=()

  while IFS=',' read -r sp r1 r2; do
    if [[ "$sp" == "$species" ]]; then
      R1_FILES+=("$r1")
      R2_FILES+=("$r2")
    fi
  done < "$TMP_FILE"

  # Concatenate files
  cat "${R1_FILES[@]}" > "$OUTPUT_DIR/${species}_R1.fastq.gz"
  cat "${R2_FILES[@]}" > "$OUTPUT_DIR/${species}_R2.fastq.gz"

  echo "→ Created $OUTPUT_DIR/${species}_R1.fastq.gz and _R2.fastq.gz"
done

# Cleanup
rm -f "$TMP_FILE"
