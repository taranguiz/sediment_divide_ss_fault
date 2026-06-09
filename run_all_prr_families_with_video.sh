#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

VIDEO_ARGS=(--with-video)

run_family() {
  local family_name="$1"
  local matrix="$2"
  local templates="$3"
  shift 3
  local labels=("$@")

  echo "Checking ${family_name} inputs..."
  python run_prr_hub_matrix.py \
    --matrix "$matrix" \
    --templates "$templates" \
    --dry-run

  for label in "${labels[@]}"; do
    echo "Starting $label with video"
    python run_prr_hub_matrix.py \
      --matrix "$matrix" \
      --templates "$templates" \
      --label "$label" \
      "${VIDEO_ARGS[@]}"
    echo "Finished $label"
  done
}

run_family \
  "Diff4" \
  "config/diffusion4_run_matrix.csv" \
  "config/diffusion4_family_templates.json" \
  Diff4_30 Diff4_20 Diff4_10 Diff4_5 Diff4_1 Diff4_05 Diff4_01

run_family \
  "TransSed" \
  "config/transport_sediment_run_matrix.csv" \
  "config/transport_sediment_family_templates.json" \
  TransSed_30 TransSed_20 TransSed_10 TransSed_5 TransSed_1 TransSed_05 TransSed_01

run_family \
  "TransSed2" \
  "config/transport_sediment_2_run_matrix.csv" \
  "config/transport_sediment_2_family_templates.json" \
  TransSed2_30 TransSed2_20 TransSed2_10 TransSed2_5 TransSed2_1 TransSed2_05 TransSed2_01

run_family \
  "TransSed3" \
  "config/transport_sediment_3_run_matrix.csv" \
  "config/transport_sediment_3_family_templates.json" \
  TransSed3_30 TransSed3_20 TransSed3_10 TransSed3_5 TransSed3_1 TransSed3_05 TransSed3_01

echo "All PRR family video runs finished."
