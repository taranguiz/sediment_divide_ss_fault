# PRR Method And Outputs

## Correct PRR Definition

Use the Duvall and Tucker profile relief ratio:

```text
PRR = R_near / R_far
```

where:

- `R_near` is maximum local relief along a profile near the fault
- `R_far` is maximum local relief along a profile farther upstream

For the current implementation, relief is measured across the full
strike-parallel profile:

```text
R = max(z profile) - min(z profile)
```

## Profile Locations

The code uses:

- fault row: `nrows / 3`
- divide/top row: `nrows - 1`
- near profile: 10% of distance from fault to divide/top boundary
- far profile: 50% of distance from fault to divide/top boundary

For the current 40 by 800 grids:

- fault row: `13`
- divide row: `39`
- near row: `16`, y = `80 m`
- far row: `26`, y = `130 m`

## Code

Main PRR helper file:

- `prr_metrics.py`

Important functions:

- `profile_row_indices`
- `compute_prr_dt`
- `compute_nae`

Postprocessing and combined plots:

- `test_prr_on_nc.py`

This scans NetCDF outputs, computes PRR, summarizes by experiment, and writes:

- `output/prr_summary/prr_dt_all_experiments.xlsx`
- `output/prr_summary/prr_dt_vs_time_all_experiments.png`
- `output/prr_summary/prr_dt_vs_nae_all_experiments.png`
- `output/prr_summary/prr_dt_vs_slip_rate_all_experiments.png`

## PRR During Runs

The dynamic loop can now sample PRR immediately after each earthquake/slip event.
This avoids needing dense NetCDF outputs.

Implemented in:

- `geomorph_dynamics_loop_trying_something.py`

Use:

```python
run_geomorf_loop(
    config,
    writer=writer,
    sample_prr_at_quakes=True,
)
```

The Duvall-Tucker 5 mm/yr test script is:

- `run_duvall_tucker_5_prr_test.py`

Smoke test:

```bash
/Users/taranguiz/opt/anaconda3/envs/landlab/bin/python run_duvall_tucker_5_prr_test.py
```

Full 100 kyr run:

```bash
RUN_FULL=1 /Users/taranguiz/opt/anaconda3/envs/landlab/bin/python run_duvall_tucker_5_prr_test.py
```

The full run output is:

- `output/Duvall_Tucker_5_prr_full/`
