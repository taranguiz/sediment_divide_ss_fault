# PRR During Run Plan

Goal: reduce storage by saving sparse NetCDF snapshots while calculating PRR
directly during the model loop after each earthquake.

Rationale:
- NetCDF snapshots are useful for inspection but are not needed for every PRR
  measurement.
- PRR after every earthquake gives a denser metric timeseries than post-hoc PRR
  from sparse NetCDF files.
- Existing 40 x 800 NetCDF snapshots are about 4 MB each, so limiting grid
  outputs to about 10 per run keeps storage low.

Implemented test path:
- `geomorph_dynamics_loop_trying_something.run_geomorf_loop(...,
  sample_prr_at_quakes=True)` samples PRR immediately after `ss_fault`.
- Outputs are written to `tabular_outputs/<model_name>/`:
  - `<model_name>_prr_at_quakes.csv`
  - `<model_name>_prr_at_quakes.xlsx`
  - `<model_name>_prr_at_quakes.png`
- `run_duvall_tucker_5_prr_test.py` runs a local Duvall_Tucker_5-style test.

Usage:
- Smoke test: `python run_duvall_tucker_5_prr_test.py`
- Full 100 kyr run: `RUN_FULL=1 python run_duvall_tucker_5_prr_test.py`

PRR definition:
- `PRR_DT = R_10 / R_50`
- `R = max(z) - min(z)` along the full strike-parallel profile.
- Near profile: 10% of fault-to-top/divide distance.
- Far profile: 50% of fault-to-top/divide distance.
