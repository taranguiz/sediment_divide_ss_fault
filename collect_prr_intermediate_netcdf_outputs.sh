#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

ARCHIVE_ROOT="download_bundles"
BUNDLE_NAME="prr_netcdf_outputs_intermediate_slip_rates"
BUNDLE_DIR="${ARCHIVE_ROOT}/${BUNDLE_NAME}"
ARCHIVE_PATH="${ARCHIVE_ROOT}/${BUNDLE_NAME}.tar.gz"

MODELS=(
  Diff4_2
  Diff4_3
  Diff4_4
  TransSed_2
  TransSed_3
  TransSed_4
  TransSed2_2
  TransSed2_3
  TransSed2_4
  TransSed3_2
  TransSed3_3
  TransSed3_4
)

rm -rf "$BUNDLE_DIR" "$ARCHIVE_PATH"
mkdir -p "$BUNDLE_DIR"

missing=0
for model in "${MODELS[@]}"; do
  src="output/${model}/netcdf_outputs"
  dest="${BUNDLE_DIR}/${model}"

  if [[ ! -d "$src" ]]; then
    echo "Missing: $src"
    missing=$((missing + 1))
    continue
  fi

  mkdir -p "$dest"
  cp -R "$src" "$dest/"
  echo "Copied $src"
done

if [[ "$missing" -ne 0 ]]; then
  echo "Stopped: $missing netcdf_outputs folder(s) were missing."
  echo "Nothing was archived. Check the missing model output folders above."
  exit 1
fi

{
  echo "PRR NetCDF outputs bundle for intermediate slip-rate runs"
  echo
  echo "Included models:"
  printf '%s\n' "${MODELS[@]}"
} > "${BUNDLE_DIR}/MANIFEST.txt"

echo "Archive size estimate before compression:"
du -sh "$BUNDLE_DIR"

tar -czf "$ARCHIVE_PATH" -C "$ARCHIVE_ROOT" "$BUNDLE_NAME"

echo
echo "Created archive:"
echo "$ARCHIVE_PATH"
echo
echo "Final archive size:"
du -sh "$ARCHIVE_PATH"
echo
echo "Download this file from the hub:"
echo "$ARCHIVE_PATH"
