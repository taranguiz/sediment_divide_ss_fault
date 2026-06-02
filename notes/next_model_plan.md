# Next Model Plan

## Existing Principle

Run new steady states only when geomorphic equilibrium parameters change.
Reuse steady states when only changing slip rate.

Parameters that require a new steady state if changed:

- `K_sed`, `K_br`
- `Vs`, `F_f`, `phi`, `H_star`
- `H0`
- `Hstar`, `V0`, `D`, `P0`, `Sc`
- `run_off`
- uplift rate, grid dimensions, resolution, or boundary conditions

Parameters that do not require a new steady state:

- `slip_rate`
- `total_model_time`
- output frequency
- PRR sampling settings

## Planned Slip-Rate Runs

These can reuse existing steady states:

| Target | Reuse Steady State |
|---|---|
| `Sed-2_5` | `final_state_Sediment_3_Duvall_Tucker.pkl` |
| `Sed-2_10` | `final_state_Sediment_3_Duvall_Tucker.pkl` |
| `Sed-2_20` | `final_state_Sediment_3_Duvall_Tucker.pkl` |
| `Sed-3_5` | `final_state_Sediment_4_Duvall_Tucker.pkl` |
| `Sed-3_20` | `final_state_Sediment_4_Duvall_Tucker.pkl` |
| `Sed-4_5` | `final_state_Sediment_5_Duvall_Tucker.pkl` |
| `Sed-4_10` | `final_state_Sediment_5_Duvall_Tucker.pkl` |
| `Sed-4_20` | `final_state_Sediment_5_Duvall_Tucker.pkl` |
| `Sed-5_05` | `final_state_Sediment_4_Duvall_Tucker_5.pkl` |
| `Sed-5_10` | `final_state_Sediment_4_Duvall_Tucker_5.pkl` |
| `Sed-5_20` | `final_state_Sediment_4_Duvall_Tucker_5.pkl` |

Optional:

| Target | Reuse Steady State |
|---|---|
| `DT_10` | `final_state_Duvall_Tucker_5.nc` |

## New Question To Hold For Future Runs

The baseline/no-sediment setup may suppress hillslope diffusion too much because
the depth-dependent diffuser needs soil thickness. The current baseline uses:

- `H0 = 0`
- `P0 = 1e-8`
- hillslope `Hstar = 1`

This may create almost no effective hillslope transport, even though the
intention was a Duvall-Tucker-like model with diffusion.

Future test idea:

- run a revised detachment-limited baseline with very small but nonzero soil
  availability, or
- use a hillslope diffuser that does not depend on soil thickness, if the goal
  is pure topographic diffusion.

Do not fold this into the planned sediment-family slip-rate runs until the
baseline design decision is made.
