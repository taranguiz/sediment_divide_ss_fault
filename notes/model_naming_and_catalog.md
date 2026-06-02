# Model Naming And Catalog

## Naming Problem

The project had three naming systems mixed together:

- original output/code names, such as `Sediment_4_Duvall_Tucker`
- manuscript/table names, such as `Sed-3_05`
- later test names, such as `Duvall_Tucker_5_prr_full`

The chosen solution is to keep original folders untouched and use manuscript
labels as the human-facing names.

## Catalog Files

See:

- `output/run_catalog/README.md`
- `output/run_catalog/current_runs.csv`
- `output/run_catalog/future_run_matrix.csv`
- `output/run_catalog/naming_scheme.md`

Alias folders were created with symlinks:

- `output/organized_runs/`
- `output/organized_steady_states/`

These are non-destructive pointers. The original folders were not renamed.

## Current Mapping

| Plot Label | Original Folder | Notes |
|---|---|---|
| `DT_05` | `output/Duvall_Tucker` | baseline, older generic steady-state file |
| `DT_5` | `output/Duvall_Tucker_5` | old 5 mm/yr baseline |
| `DT_5_new` | `output/Duvall_Tucker_5_prr_full` | new lower-domain-moving 100 kyr run |
| `Sed-1_05` | `output/Sediment_2_Duvall_Tucker` | fluvial sediment setup |
| `Sed-2_05` | `output/Sediment_3_Duvall_Tucker` | lower runoff |
| `Sed-3_05` | `output/Sediment_4_Duvall_Tucker` | active hillslope setup |
| `Sed-3_10` | `output/Sediment_4_Duvall_Tucker_20_` | active hillslope, 10 mm/yr |
| `Sed-4_05` | `output/Sediment_5_Duvall_Tucker` | higher soil production |
| `Sed-5_5` | `output/Sediment_4_Duvall_Tucker_5` | no-fines, 5 mm/yr |
| `INTERMEDIATE` | `output/Sediment_Duvall_Tucker` | skipped in main plot |

## Recommended Future Names

Use manuscript labels in plots:

```text
Sed-2_10
Sed-3_20
Sed-4_5
Sed-5_05
```

Use machine-safe folder/model names in code:

```text
Sed_2_10
Sed_3_20
Sed_4_5
Sed_5_05
```

Use `05` for `0.5 mm/yr`; use plain integers for `5`, `10`, and `20 mm/yr`.
