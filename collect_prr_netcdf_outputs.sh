#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

ARCHIVE_ROOT="download_bundles"
BUNDLE_NAME="prr_netcdf_outputs_all_families"
BUNDLE_DIR="${ARCHIVE_ROOT}/${BUNDLE_NAME}"
ARCHIVE_PATH="${ARCHIVE_ROOT}/${BUNDLE_NAME}.tar.gz"

MODELS=(
  Diff4_01
  Diff4_05
  Diff4_1
  Diff4_5
  Diff4_10
  Diff4_20
  Diff4_30
  TransSed_01
  TransSed_05
  TransSed_1
  TransSed_5
  TransSed_10
  TransSed_20
  TransSed_30
  TransSed2_01
  TransSed2_05
  TransSed2_1
  TransSed2_5
  TransSed2_10
  TransSed2_20
  TransSed2_30
  TransSed3_01
  TransSed3_05
  TransSed3_1
  TransSed3_5
  TransSed3_10
  TransSed3_20
  TransSed3_30
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
  echo "PRR NetCDF outputs bundle"
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
