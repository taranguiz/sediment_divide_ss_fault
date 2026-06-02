# Sediment Divide Strike-Slip Fault Models

This repository contains Landlab scripts for synthetic divide/ridge mobility
experiments around a right-lateral strike-slip fault. The current workflow
focuses on profile relief ratio (PRR) measured after earthquake/slip events and
compared across model families and slip rates.

## Current Analysis Goal

The near-term goal is to run slip-rate variants for the sediment model families
and compare:

- PRR vs slip rate
- PRR vs `Nae`
- PRR time series after earthquake events

Runs in `config/prr_hub_run_matrix.csv` use the same total_slip,
currently `1000 m`. The runner computes `total_model_time` from slip rate, so
slow faults run longer than fast faults.

The current PRR definition is:

```text
PRR = R_near / R_far
```

where `R_near` is full strike-parallel profile relief at 10% of the
fault-to-divide distance, and `R_far` is full strike-parallel profile relief at
50% of that distance.

## Repository Layout

- `geomorph_dynamics_loop_trying_something.py` - dynamic faulting loop.
- `updated_steady.py` - steady-state builder.
- `ss_fault_function.py` - strike-slip grid motion.
- `prr_metrics.py` - PRR and `Nae` helper functions.
- `run_prr_hub_matrix.py` - hub-ready runner for planned PRR experiments.
- `config/prr_hub_run_matrix.csv` - planned slip-rate matrix.
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

By default, the hub runner saves sparse NetCDF snapshots and PRR/event tables,
but skips topography PNG frames to reduce disk use. To save PNG frames too:

```bash
python run_prr_hub_matrix.py --label Sed-3_20 --with-topo-plots
```

The matrix stores `total_slip`, not fixed run duration. For the current
`total_slip = 1000 m`, run durations are:

```text
0.5 mm/yr -> 2,000,000 yr
5 mm/yr   ->   200,000 yr
10 mm/yr  ->   100,000 yr
20 mm/yr  ->    50,000 yr
```

## Required Steady States For Planned Runs

Upload these files to `output/steady_state_files/`:

```text
final_state_Sediment_3_Duvall_Tucker.pkl
final_state_Sediment_4_Duvall_Tucker.pkl
final_state_Sediment_5_Duvall_Tucker.pkl
final_state_Sediment_4_Duvall_Tucker_5.pkl
```

Optional detachment-limited comparison:

```text
final_state_Duvall_Tucker_5.nc
```

## Notes

Start with `notes/README.md` when returning to the project. It points to the
recommended reading order for the current modeling decisions and naming scheme.
