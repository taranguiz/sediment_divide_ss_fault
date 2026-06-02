# Session Handoff: 2026-06-02

## Current Project Goal

The project is modeling divide/ridge mobility around strike-slip faults using
Landlab. The current analysis compares profile relief ratio (PRR) across
different model families and slip rates, with special attention to whether
sediment/fluvial/hillslope parameters change ridge mobility.

The immediate next goal is to run additional slip-rate variants to populate the
combined PRR vs `Nae` and PRR vs slip-rate plots.

## Major Decisions So Far

- Use the Duvall and Tucker PRR definition:
  `PRR = near-fault relief / farther-upstream relief`.
- Use full strike-parallel profile relief, not a moving-window/local relief
  method, for the current plots.
- The near-fault profile is at 10% of the distance from the fault to the top
  boundary.
- The far-upstream profile is at 50% of the distance from the fault to the top
  boundary.
- Mark sediment-effect models based on `K_sed != K_br`.
- Skip `Sediment_Duvall_Tucker` in the main plot because it was an intermediate
  run and not part of the manuscript table.
- Keep old output folders intact; use catalogs and symlinks to organize names.
- For future slip-rate-only variants, reuse the matching steady-state grid.
- If geomorphic equilibrium parameters change, run a new steady state first.

## Important Modeling Issue To Remember

The Duvall-Tucker-style baseline tried to suppress sediment by using:

- `H_star = 1e-8` in SPACE
- `Vs = 0`
- `F_f = 1`
- `H0 = 0`
- very low `P0 = 1e-8`

This makes sense for making the fluvial/channel part close to bare-bedrock or
detachment-limited behavior. However, the hillslope component is
depth-dependent. With `H0 = 0` and very low `P0`, soil thickness may remain near
zero, which can suppress the effective hillslope diffusion. This should be held
as a possible revision to the future baseline model plan: a small amount of soil
may be needed if we want hillslope diffusion to actually operate.

## Current Clean Outputs

Combined PRR outputs:

- `output/prr_summary/prr_dt_all_experiments.xlsx`
- `output/prr_summary/prr_dt_vs_time_all_experiments.png`
- `output/prr_summary/prr_dt_vs_nae_all_experiments.png`
- `output/prr_summary/prr_dt_vs_slip_rate_all_experiments.png`

New 100 kyr Duvall-Tucker 5 mm/yr run with updated fault motion:

- output folder: `output/Duvall_Tucker_5_prr_full/`
- video: `output/Duvall_Tucker_5_prr_full/Duvall_Tucker_5_prr_full_evolution.mp4`
- PRR after every earthquake:
  `output/Duvall_Tucker_5_prr_full/tabular_outputs/Duvall_Tucker_5_prr_full/Duvall_Tucker_5_prr_full_prr_at_quakes.csv`

## Current Key Result

`DT_5_new` was added to the combined PRR plot as the new run with the lower
domain moving:

- slip rate: `5 mm/yr`
- `Nae`: `13.354701`
- PRR mean: `0.868996`
- PRR std: `0.089022`
- measurements: `100`, sampled after earthquakes

For comparison, the older `DT_5` sparse-NetCDF result was:

- PRR mean: `0.919365`
- PRR std: `0.066874`
- measurements: `11`, sampled from NetCDF snapshots

## Resource Estimate From Local Test

The 100 kyr sparse-output test used:

- about `47 MB` total disk
- `10` NetCDF snapshots
- one PRR measurement after each earthquake

Recommended hub request:

- memory per run: `1-2 GB`
- recommended allocation: `2 GB/run`
- on a 4 GB hub: run one model at a time
- disk estimate: about `100-150 MB/run`; request several GB for safety
