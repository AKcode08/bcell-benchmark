#!/bin/bash
#SBATCH --job-name=graphbepi_batch
#SBATCH --output=logs/graphbepi_batch_%j.out
#SBATCH --error=logs/graphbepi_batch_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=96G
#SBATCH --cpus-per-task=8

module purge
module load bluebear
module load bear-apps/2024a/live
module load Miniforge3/25.3.0-3

source $(conda info --base)/etc/profile.d/conda.sh
conda activate graphbepi

PROJECT_DIR="/rds/projects/e/elhamsak-epitope-dev/bcell_benchmark"
GRAPHBEPI_DIR="$PROJECT_DIR/models/GraphBepi"

# Give pathogens like: ESKAPE,GAST,INTRA,RESP
PATHOGENS=${PATHOGENS:-"ESKAPE"}
TYPES=${TYPES:-"linear,conformational"}

mkdir -p "$PROJECT_DIR/logs/graphbepi"

echo "Using Python:"
which python

echo "Checking GraphBepi environment:"
python -c "import esm, torch, pytorch_lightning; print('Environment OK')"

cd "$GRAPHBEPI_DIR" || exit 1

IFS=',' read -ra PATHOGEN_LIST <<< "$PATHOGENS"
IFS=',' read -ra TYPE_LIST <<< "$TYPES"

for pathogen in "${PATHOGEN_LIST[@]}"; do
  for epitype in "${TYPE_LIST[@]}"; do

    INPUT_DIR="$PROJECT_DIR/structures/$pathogen/$epitype"
    OUTPUT_DIR="$PROJECT_DIR/outputs/graphbepi/$pathogen/$epitype"
    TMP_DIR="$OUTPUT_DIR/tmp"

    mkdir -p "$OUTPUT_DIR"
    mkdir -p "$TMP_DIR"

    echo "======================================"
    echo "Pathogen: $pathogen"
    echo "Type: $epitype"
    echo "Input: $INPUT_DIR"
    echo "Output: $OUTPUT_DIR"
    echo "======================================"

    for pdb in "$INPUT_DIR"/*.pdb; do

      [ -e "$pdb" ] || continue

      base=$(basename "$pdb" .pdb)
      final_csv="$OUTPUT_DIR/${base}_graphbepi.csv"
      antigen_tmp="$TMP_DIR/$base"

      if [ -f "$final_csv" ]; then
        echo "Skipping $base — output already exists."
        continue
      fi

      echo "Running GraphBepi on $base"

      rm -rf "$antigen_tmp"
      mkdir -p "$antigen_tmp"

      python test.py \
        -i "$pdb" \
        -p \
        --gpu -1 \
        -o "$antigen_tmp"

      csv_file=$(find "$antigen_tmp" -name "*.csv" | head -n 1)

      if [ -f "$csv_file" ]; then
        mv "$csv_file" "$final_csv"
        echo "Saved: $final_csv"
      else
        echo "WARNING: No CSV produced for $base"
      fi

      rm -rf "$antigen_tmp"

    done

    rmdir "$TMP_DIR" 2>/dev/null || true

  done
done

echo "GraphBepi batch finished."
