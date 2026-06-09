#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

MATRIX_2="config/transport_sediment_2_run_matrix.csv"
TEMPLATES_2="config/transport_sediment_2_family_templates.json"
MATRIX_3="config/transport_sediment_3_run_matrix.csv"
TEMPLATES_3="config/transport_sediment_3_family_templates.json"
VIDEO_ARGS=(--with-video)

echo "Checking TransSed2 inputs..."
python run_prr_hub_matrix.py \
  --matrix "$MATRIX_2" \
  --templates "$TEMPLATES_2" \
  --dry-run

echo "Checking TransSed3 inputs..."
python run_prr_hub_matrix.py \
  --matrix "$MATRIX_3" \
  --templates "$TEMPLATES_3" \
  --dry-run

for label in TransSed2_30 TransSed2_20 TransSed2_10 TransSed2_5 TransSed2_1 TransSed2_05 TransSed2_01; do
  echo "Starting $label"
  python run_prr_hub_matrix.py \
    --matrix "$MATRIX_2" \
    --templates "$TEMPLATES_2" \
    --label "$label" \
    "${VIDEO_ARGS[@]}"
  echo "Finished $label"
done

for label in TransSed3_30 TransSed3_20 TransSed3_10 TransSed3_5 TransSed3_1 TransSed3_05 TransSed3_01; do
  echo "Starting $label"
  python run_prr_hub_matrix.py \
    --matrix "$MATRIX_3" \
    --templates "$TEMPLATES_3" \
    --label "$label" \
    "${VIDEO_ARGS[@]}"
  echo "Finished $label"
done

echo "All requested TransSed2 and TransSed3 video runs finished."
