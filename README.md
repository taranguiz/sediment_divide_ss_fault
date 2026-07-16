# Sediment Divide Strike-Slip Fault Models

This repository contains Landlab scripts for synthetic divide/ridge mobility
experiments around a right-lateral strike-slip fault. The current workflow
focuses on extracting profile relief ratio (PRR) (Duvall and Tucker, 2015), measured after earthquake/slip events and
compared across model families and slip rates.

## Repository Layout

- `geomorph_dynamics_loop_trying_something.py` - dynamic faulting loop.
- `updated_steady.py` - steady-state builder.
- `ss_fault_function.py` - strike-slip grid motion.
- `prr_metrics.py` - PRR and `Nae` helper functions.
- `run_prr_hub_matrix.py` - hub-ready runner for planned PRR experiments.
- `config/prr_hub_run_matrix.csv` - planned slip-rate matrix.
- `config/prr_hub_run_matrix_all.csv` - comprehensive official rerun matrix.
- `config/prr_family_templates.json` - committed family parameter templates.
- `notes/` - project notes, handoff context, naming rules, and PRR method notes.

Generated outputs are intentionally ignored by Git:

- `output/`
- `output_topo/`
- `__pycache__/`
- `.DS_Store`

## Hub Workflow

Do not commit NetCDF, pickle, or generated output files. Clone this repository
on the hub, then upload/copy required steady-state files into:

```text
output/steady_state_files/
```

Check that the selected runs can see their steady-state files:

```bash
python run_prr_hub_matrix.py --dry-run
```

Run one model:

```bash
python run_prr_hub_matrix.py --label Sed-3_20
```

Run all enabled planned models:

```bash
python run_prr_hub_matrix.py
```

Run the comprehensive official rerun matrix, including all families

By default, the hub runner saves sparse NetCDF snapshots and PRR/event tables,
but skips topography PNG frames to reduce disk use. To save PNG frames too:

```bash
python run_prr_hub_matrix.py --label Sed-3_20 --with-topo-plots
```

To save an MP4 evolution video:

```bash
python run_prr_hub_matrix.py --label Sed-3_20 --with-video
```

`--with-video` automatically saves the sparse topography PNG frames used to
build the movie and deletes those frame PNGs after they are added to the MP4.
Use `--keep-video-frames` if you also want to keep the PNGs.

The matrix stores `total_slip`, not fixed run duration. For the current
`total_slip = 1000 m`, run durations are:

```text
0.5 mm/yr -> 2,000,000 yr
5 mm/yr   ->   200,000 yr
10 mm/yr  ->   100,000 yr
20 mm/yr  ->    50,000 yr
```

## Required Steady States For Planned Runs

Upload these files to `output/steady_state_files/`
