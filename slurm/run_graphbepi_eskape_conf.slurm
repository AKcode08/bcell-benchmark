#!/bin/bash
#SBATCH --job-name=gb_esk_conf
#SBATCH --output=logs/gb_esk_conf.out
#SBATCH --error=logs/gb_esk_conf.err
#SBATCH --time=12:00:00
#SBATCH --mem=96G
#SBATCH --cpus-per-task=8

module purge
module load bluebear
module load bear-apps/2024a/live
module load Miniforge3/25.3.0-3

source $(conda info --base)/etc/profile.d/conda.sh
conda activate graphbepi

echo "Python being used:"
which python

echo "Checking esm install:"
python -c "import esm; print('esm works')"

PROJECT_DIR="/rds/projects/e/elhamsak-epitope-dev/bcell_benchmark"
GRAPHBEPI_DIR="$PROJECT_DIR/models/GraphBepi"
INPUT_DIR="$PROJECT_DIR/structures/ESKAPE/conformational"
OUTPUT_DIR="$PROJECT_DIR/outputs/graphbepi/ESKAPE/conformational"

mkdir -p "$PROJECT_DIR/logs"
mkdir -p "$OUTPUT_DIR"

cd "$GRAPHBEPI_DIR" || exit 1

for pdb in "$INPUT_DIR"/*.pdb; do

    base=$(basename "$pdb" .pdb)

    echo "======================================"
    echo "Running GraphBepi on: $base"
    echo "======================================"

    python test.py \
        -i "$pdb" \
        -p \
        --gpu -1 \
        -o "$OUTPUT_DIR/$base"

done

echo "Batch finished."