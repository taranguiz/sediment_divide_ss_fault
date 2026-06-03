# OpenEarthscape Hub Setup

Hub URL:

```text
https://frontier.openearthscape.org/user/taranguiz/lab
```

## GitHub

The local repository remote is:

```text
https://github.com/taranguiz/sediment_divide_ss_fault.git
```

On the hub, clone or pull this repository, then work from the repository root.

## Data Policy

Keep generated data out of Git. This includes:

- NetCDF snapshots
- steady-state pickle/NetCDF files
- output plots and videos
- PRR summary workbooks generated from runs

These files belong in hub storage or local storage, not in Git history.

## Required Steady-State Uploads

Create this folder on the hub:

```bash
mkdir -p output/steady_state_files
```

Upload/copy these files into it:

```text
final_state_Sediment_2_Duvall_Tucker.pkl
final_state_Sediment_3_Duvall_Tucker.pkl
final_state_Sediment_4_Duvall_Tucker.pkl
final_state_Sediment_5_Duvall_Tucker.pkl
final_state_Sediment_4_Duvall_Tucker_5.pkl
final_state_Duvall_Tucker_5.nc
```

The comprehensive official rerun matrix needs all six files above. Together
they are about 38 MB locally.

## Check Before Running

From the repository root:

```bash
python run_prr_hub_matrix.py --dry-run
```

This checks that each selected run has its steady-state file available.

## Running Models

Run one model at a time on a small hub allocation:

```bash
python run_prr_hub_matrix.py --label Sed-3_20
```

Run all enabled sediment variants:

```bash
python run_prr_hub_matrix.py
```

Run the comprehensive official rerun matrix, including DT and Sed-1 through
Sed-5:

```bash
python run_prr_hub_matrix.py --matrix config/prr_hub_run_matrix_all.csv
```

## Outputs

Each run writes to:

```text
output/<model_name>/
```

Important PRR files:

```text
output/<model_name>/tabular_outputs/<model_name>/<model_name>_prr_at_quakes.csv
output/<model_name>/tabular_outputs/<model_name>/<model_name>_prr_at_quakes.xlsx
```

Sparse NetCDF snapshots:

```text
output/<model_name>/netcdf_outputs/<model_name>/
```

The hub runner keeps topography PNG frames off by default. This helps keep the
run close to the goal of about 10 NetCDF snapshots plus tabular PRR outputs.
To save an MP4 evolution video for a specific run, add `--with-video`:

```bash
python run_prr_hub_matrix.py --matrix config/prr_hub_run_matrix_all.csv --label Sed-3_20 --with-video
```

By default, `--with-video` deletes each temporary PNG frame after adding it to
the MP4. Use `--keep-video-frames` if you also want the frame PNGs.

The comprehensive matrix writes output folders prefixed with `PRR1000_`, such
as `output/PRR1000_Sed_3_20/`. This keeps official 1000 m total-slip reruns
separate from earlier test outputs.

The run matrix uses the same `total_slip` for every slip rate. The
current `total_slip` is `1000 m`, and `run_prr_hub_matrix.py` computes run duration
from:

```text
total_model_time = total_slip / (slip_rate_mm_yr / 1000)
```

For `total_slip = 1000 m`:

```text
0.5 mm/yr -> 2,000,000 yr
5 mm/yr   ->   200,000 yr
10 mm/yr  ->   100,000 yr
20 mm/yr  ->    50,000 yr
```

## Planned Naming

Use manuscript labels in plots:

```text
Sed-3_20
```

Use machine-safe model/output folder names:

```text
Sed_3_20
```

Slip-rate encoding:

```text
05 = 0.5 mm/yr
5 = 5 mm/yr
10 = 10 mm/yr
20 = 20 mm/yr
```
