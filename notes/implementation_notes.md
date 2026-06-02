# Implementation Notes

## Compatibility Helpers

Added:

- `numpy_compat.py`
- `landlab_compat.py`

These address compatibility issues with the local Landlab/NumPy setup.

## Fault Motion

`ss_fault_function.py` was modified so the lower/bottom side of the domain moves
instead of shifting the top side. This was important for the new
`Duvall_Tucker_5_prr_full` comparison run.

## Dynamic Loop Changes

File:

- `geomorph_dynamics_loop_trying_something.py`

Changes:

- `read_grid(config)` now supports `initial_state_model_name`.
- This allows a new output model to reuse another model's steady-state grid.
- Added optional PRR sampling after each earthquake:
  `sample_prr_at_quakes=True`.
- Added sparse-output-friendly behavior.
- Added optional `writer=None`, `interactive_plots`, and `save_outputs`
  controls.
- PRR-at-quakes outputs are saved as CSV, XLSX, and PNG.

## Steady-State Building

File:

- `updated_steady.py`

Current behavior:

- builds steady state from scratch using current config parameters
- writes model-specific steady file:
  `output/steady_state_files/final_state_<model_name>.pkl`

Important rule:

- new geomorphic parameters require a new steady state
- slip-rate-only variants can reuse an existing steady state

## Combined Plot Changes

File:

- `test_prr_on_nc.py`

Changes:

- Uses corrected `PRR_DT = R_10 / R_50`.
- Adds table/manuscript labels.
- Skips `Sediment_Duvall_Tucker`.
- Marks sediment-effect runs using `K_sed != K_br`.
- Includes additional precomputed PRR time series:
  `Duvall_Tucker_5_prr_full` as `DT_5_new`.

## Known Cleanup Issues

- The dynamic loop is still very verbose during runs.
- There is an old mean-timeseries length mismatch warning. PRR outputs are not
  affected.
- `DT_05` is older and uses generic `final_state.pkl`; future baseline runs
  should use named steady-state files.
- `Sed-3_10` has its own steady-state file, but it is identical to `Sed-3_05`;
  this is redundant, not wrong.
