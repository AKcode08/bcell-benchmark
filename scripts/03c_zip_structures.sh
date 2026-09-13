#!/usr/bin/env bash
set -euo pipefail

BASE="structures"

for ORG in "$BASE"/*; do
  for TYPE in linear conformational; do

    DIR="$ORG/$TYPE"
    [[ -d "$DIR" ]] || continue

    echo "Processing: $DIR"

    # Get all pdb files (portable)
    FILES=$(find "$DIR" -maxdepth 1 -name "*.pdb" | sort)

    TOTAL=$(echo "$FILES" | grep -c . || true)
    echo "Total PDBs: $TOTAL"

    if [[ "$TOTAL" -eq 0 ]]; then
      continue
    fi

    # Convert to array safely
    IFS=$'\n' read -d '' -r -a FILE_ARRAY <<< "$FILES" || true

    # If ≤100 → single zip
    if [[ "$TOTAL" -le 100 ]]; then
      ZIPFILE="$DIR/$(basename "$ORG")_${TYPE}.zip"
      zip -j "$ZIPFILE" "${FILE_ARRAY[@]}"
      echo "Created: $ZIPFILE"

    else
      PART=1
      COUNT=0
      CHUNK=()

      for f in "${FILE_ARRAY[@]}"; do
        CHUNK+=("$f")
        ((COUNT++))

        if [[ "$COUNT" -eq 100 ]]; then
          ZIPFILE="$DIR/$(basename "$ORG")_${TYPE}_part${PART}.zip"
          zip -j "$ZIPFILE" "${CHUNK[@]}"
          echo "Created: $ZIPFILE"

          CHUNK=()
          COUNT=0
          ((PART++))
        fi
      done

      # Remaining files
      if [[ "$COUNT" -gt 0 ]]; then
        ZIPFILE="$DIR/$(basename "$ORG")_${TYPE}_part${PART}.zip"
        zip -j "$ZIPFILE" "${CHUNK[@]}"
        echo "Created: $ZIPFILE"
      fi
    fi

    echo "--------------------------------------"

  done
done