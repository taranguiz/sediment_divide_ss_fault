#!/usr/bin/env python3
"""Compute PRR through time from NetCDF experiment outputs.

PRR is calculated from swath-averaged strike-parallel profiles:

- near-fault swath: just upstream of the fault to 10% of fault-to-divide distance
- far swath: 40% to 50% of fault-to-divide distance
- PRR = relief(near swath-averaged profile) / relief(far swath-averaged profile)

For long domains, this script uses the full x-domain. For other domains, it
uses a centered 1000 m segment when the domain is longer than 1000 m.
"""
from __future__ import annotations

import glob
import argparse
import json
import os
import re
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import numpy_compat  # noqa: F401
import pandas as pd
from landlab.io.netcdf import read_netcdf

from prr_metrics import (
    compute_nae,
    compute_prr_swath_profile,
    far_fault_swath_row_indices,
    near_fault_swath_row_indices,
    x_slice_from_segment,
)


DEFAULT_INPUT_ROOTS = ",".join(
    [
        "/Volumes/T7Shield/Ch4_outputs",
        "/Users/taranguiz/Research/sediment_divide_ss_fault/output/Duvall_Tucker",
        "/Users/taranguiz/Research/sediment_divide_ss_fault/output/Duvall_Tucker_5/netcdf_outputs/Duvall_Tucker_5",
    ]
)
DEFAULT_LOCAL_CONFIG_ROOT = "/Users/taranguiz/Research/sediment_divide_ss_fault/output"
DEFAULT_RESULTS_DIR = "/Users/taranguiz/Research/sediment_divide_ss_fault/output/prr_summary"
DEFAULT_EXCLUDE_EXPERIMENTS = {"Sediment_Duvall_Tucker"}
EXPERIMENT_LABELS = {
    "Duvall_Tucker": "DT_05",
    "Duvall_Tucker_5": "DT_5",
    "Sediment_2_Duvall_Tucker": "Sed-1_05",
    "Sediment_3_Duvall_Tucker": "Sed-2_05",
    "Sediment_4_Duvall_Tucker": "Sed-3_05",
    "Sediment_4_Duvall_Tucker_20_": "Sed-3_10",
    "Sediment_5_Duvall_Tucker": "Sed-4_05",
    "Sediment_4_Duvall_Tucker_5": "Sed-5_5",
    "S01_Cont_big_domain": "S01_big_0.1",
    "S05_Cont_big_domain": "S05_big_0.5",
}


def _prompt_path(prompt: str, default_path: str) -> str:
    raw = input(f"{prompt} [{default_path}]: ").strip()
    return raw or default_path


def _path_arg_or_prompt(value: str | None, prompt: str, default_path: str) -> str:
    if value:
        return value
    if sys.stdin.isatty():
        return _prompt_path(prompt, default_path)
    return default_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compute PRR through time from NetCDF experiment outputs using the "
            "swath-profile PRR method."
        )
    )
    parser.add_argument(
        "--input-roots",
        help="Comma-separated experiment folders or parent folders containing NetCDF outputs.",
    )
    parser.add_argument(
        "--config-root",
        help="Local root used as a fallback when config.json is not inside an experiment folder.",
    )
    parser.add_argument(
        "--results-dir",
        help="Directory where PRR summary tables and plots will be saved.",
    )
    return parser.parse_args()


def _extract_time_from_filename(path: str) -> int:
    match = re.search(r"(\d+)\.nc$", os.path.basename(path))
    return int(match.group(1)) if match else -1


def _load_experiment_config(exp_dir: str, local_config_root: str) -> dict | None:
    exp_name = os.path.basename(exp_dir.rstrip("/"))
    candidate_paths = [
        os.path.join(exp_dir, "config.json"),
        os.path.join(local_config_root, exp_name, "config.json"),
    ]
    for cfg_path in candidate_paths:
        if os.path.exists(cfg_path):
            with open(cfg_path, "r") as f:
                cfg = json.load(f)
            cfg["_config_path"] = cfg_path
            return cfg
    return None


def _experiment_name_from_config(cfg: dict, exp_dir: str) -> str:
    return str(cfg.get("model_name") or os.path.basename(exp_dir.rstrip("/")))


def _slip_rate_from_config(cfg: dict) -> float:
    if cfg.get("slip_rate") is not None:
        return float(cfg["slip_rate"])
    total_slip = cfg.get("total_slip")
    total_model_time = cfg.get("total_model_time")
    if total_slip is None or total_model_time in (None, 0):
        raise KeyError(
            "Config has no slip_rate and cannot infer it from total_slip / total_model_time"
        )
    return (float(total_slip) / float(total_model_time)) * 1000.0


def _netcdf_files(exp_dir: str) -> list[str]:
    return sorted(
        path for path in glob.glob(os.path.join(exp_dir, "*.nc"))
        if not os.path.basename(path).startswith("._")
    )


def _split_input_roots(raw: str) -> list[str]:
    return [part.strip() for part in raw.split(",") if part.strip()]


def _discover_experiment_dirs(input_roots: list[str]) -> list[str]:
    exp_dirs = []
    seen = set()
    for root in input_roots:
        if not os.path.isdir(root):
            print(f"Skipping missing input root: {root}")
            continue
        candidates = [root] if _netcdf_files(root) else [
            d for d in glob.glob(os.path.join(root, "*")) if os.path.isdir(d)
        ]
        for candidate in candidates:
            resolved = os.path.realpath(candidate)
            if resolved in seen:
                continue
            if _netcdf_files(candidate):
                exp_dirs.append(candidate)
                seen.add(resolved)
    return sorted(exp_dirs)


def _x_slice_for_experiment(model_name: str, ncols: int, dxy: float) -> tuple[slice, str]:
    if "big_domain" in model_name:
        return x_slice_from_segment(ncols, dxy), "full_domain"
    return x_slice_from_segment(ncols, dxy, segment_length_m=1000.0), "centered_1000m"


def _annotate_points(ax, df: pd.DataFrame, x_col: str, y_col: str) -> None:
    for _, row in df.iterrows():
        ax.annotate(
            row["plot_label"],
            (row[x_col], row[y_col]),
            xytext=(5, 5),
            textcoords="offset points",
            fontsize=8,
        )


def _plot_summary(
    df_summary: pd.DataFrame,
    x_col: str,
    x_label: str,
    output_path: str,
) -> None:
    fig, ax = plt.subplots(figsize=(8, 6))
    for has_sediment_effect, grp in df_summary.groupby("has_sediment_effect", sort=True):
        ax.errorbar(
            grp[x_col],
            grp["PRR_mean_time"],
            yerr=grp["PRR_std_time"],
            fmt="s" if has_sediment_effect else "o",
            capsize=4,
            linestyle="none",
            label="sediment effect" if has_sediment_effect else "uniform erodibility",
        )
    _annotate_points(ax, df_summary, x_col, "PRR_mean_time")
    ax.axhline(1.0, color="gray", ls="--", lw=0.8)
    ax.set_xlabel(x_label)
    ax.set_ylabel("PRR")
    ax.set_title(f"PRR vs {x_label}")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    fig.savefig(output_path, dpi=220, facecolor="white")
    plt.close(fig)
    print(f"Saved plot: {output_path}")


def main() -> None:
    args = parse_args()
    input_root_raw = _path_arg_or_prompt(
        args.input_roots,
        "Enter experiment output directories (comma-separated)",
        DEFAULT_INPUT_ROOTS,
    )
    input_roots = _split_input_roots(input_root_raw)
    local_config_root = _path_arg_or_prompt(
        args.config_root,
        "Enter local config root directory",
        DEFAULT_LOCAL_CONFIG_ROOT,
    )
    results_dir = _path_arg_or_prompt(
        args.results_dir,
        "Enter directory to save PRR results",
        DEFAULT_RESULTS_DIR,
    )

    exp_dirs = _discover_experiment_dirs(input_roots)
    if not exp_dirs:
        sys.exit(f"No experiment folders with .nc files found under: {input_roots}")

    print(f"Found {len(exp_dirs)} experiment folders.")
    all_rows = []
    summary_rows = []

    for exp_dir in exp_dirs:
        exp_name = os.path.basename(exp_dir.rstrip("/"))
        cfg = _load_experiment_config(exp_dir, local_config_root)
        if cfg is None:
            print(f"Skipping {exp_name}: config.json not found.")
            continue

        model_name = _experiment_name_from_config(cfg, exp_dir)
        if model_name in DEFAULT_EXCLUDE_EXPERIMENTS:
            print(f"Skipping {model_name}: excluded from default PRR summary.")
            continue

        nc_files = _netcdf_files(exp_dir)
        if not nc_files:
            continue

        slip_rate = _slip_rate_from_config(cfg)
        k_sed = float(cfg["K_sed"])
        k_br = float(cfg["K_br"])
        d_coef = float(cfg["D"])
        total_steady_time = int(cfg.get("total_steady_time", 0))
        has_sediment_effect = not np.isclose(k_sed, k_br)
        plot_label = EXPERIMENT_LABELS.get(model_name, model_name)

        first_grid = read_netcdf(nc_files[0])
        nrows = int(first_grid.number_of_node_rows)
        ncols = int(first_grid.number_of_node_columns)
        dxy = float(first_grid.dx)
        fault_row = nrows // 3
        divide_row = nrows - 1
        row_near_start, row_near_end = near_fault_swath_row_indices(fault_row, divide_row)
        row_far_start, row_far_end = far_fault_swath_row_indices(fault_row, divide_row)
        x_slice, prr_x_method = _x_slice_for_experiment(model_name, ncols, dxy)
        nae = compute_nae(slip_rate, k_br, d_coef)

        print(
            f"\n{exp_name}: {len(nc_files)} snapshots | "
            f"near rows {row_near_start}-{row_near_end} | "
            f"far rows {row_far_start}-{row_far_end} | "
            f"x method={prr_x_method} | grid={nrows}x{ncols} | "
            f"Nae={nae:.4f} | config={cfg['_config_path']}"
        )

        exp_rows = []
        for i, nc_path in enumerate(nc_files):
            calendar_yr = _extract_time_from_filename(nc_path)
            model_yr = calendar_yr - total_steady_time if calendar_yr >= 0 else np.nan
            mg = read_netcdf(nc_path)
            file_nrows = int(mg.number_of_node_rows)
            file_ncols = int(mg.number_of_node_columns)
            if (file_nrows, file_ncols) != (nrows, ncols):
                raise ValueError(
                    f"{nc_path} has grid {file_nrows}x{file_ncols}, expected {nrows}x{ncols}"
                )

            z_2d = mg.at_node["topographic__elevation"].reshape((file_nrows, file_ncols))
            prr = compute_prr_swath_profile(
                z_2d,
                row_near_start,
                row_near_end,
                row_far_start,
                row_far_end,
                x_slice=x_slice,
            )
            row = {
                "experiment": exp_name,
                "model_name": model_name,
                "plot_label": plot_label,
                "snapshot": i,
                "calendar_yr": calendar_yr,
                "model_yr": model_yr,
                "Nae": nae,
                "slip_rate_mm_yr": slip_rate,
                "K_sed": k_sed,
                "K_br": k_br,
                "has_sediment_effect": has_sediment_effect,
                "D": d_coef,
                "dxy": dxy,
                "fault_row": fault_row,
                "divide_row": divide_row,
                "row_near_swath_start": row_near_start,
                "row_near_swath_end": row_near_end,
                "row_far_swath_start": row_far_start,
                "row_far_swath_end": row_far_end,
                "fault_y": float(first_grid.node_y[fault_row * ncols]),
                "divide_y": float(first_grid.node_y[divide_row * ncols]),
                "y_near_swath_start": float(first_grid.node_y[row_near_start * ncols]),
                "y_near_swath_end": float(first_grid.node_y[row_near_end * ncols]),
                "y_far_swath_start": float(first_grid.node_y[row_far_start * ncols]),
                "y_far_swath_end": float(first_grid.node_y[row_far_end * ncols]),
                "prr_x_method": prr_x_method,
                "x_col_start": prr["x_col_start"],
                "x_col_end_exclusive": prr["x_col_end_exclusive"],
                "n_x_cols": prr["n_x_cols"],
                "PRR": prr["PRR_swath_profile"],
                "R_near_swath_profile": prr["R_near_swath_profile"],
                "R_far_swath_profile": prr["R_far_swath_profile"],
                "measurement_source": "netcdf_snapshots",
            }
            exp_rows.append(row)
            all_rows.append(row)

        exp_df = pd.DataFrame(exp_rows)
        summary_rows.append(
            {
                "experiment": exp_name,
                "model_name": model_name,
                "plot_label": plot_label,
                "Nae": nae,
                "slip_rate_mm_yr": slip_rate,
                "K_sed": k_sed,
                "K_br": k_br,
                "has_sediment_effect": has_sediment_effect,
                "D": d_coef,
                "dxy": dxy,
                "fault_row": fault_row,
                "divide_row": divide_row,
                "row_near_swath_start": row_near_start,
                "row_near_swath_end": row_near_end,
                "row_far_swath_start": row_far_start,
                "row_far_swath_end": row_far_end,
                "prr_x_method": prr_x_method,
                "n_snapshots": len(exp_df),
                "measurement_source": "netcdf_snapshots",
                "PRR_mean_time": float(exp_df["PRR"].mean()),
                "PRR_std_time": float(exp_df["PRR"].std()),
                "PRR_min_time": float(exp_df["PRR"].min()),
                "PRR_max_time": float(exp_df["PRR"].max()),
                "R_near_swath_profile_mean": float(exp_df["R_near_swath_profile"].mean()),
                "R_far_swath_profile_mean": float(exp_df["R_far_swath_profile"].mean()),
            }
        )

    if not all_rows:
        sys.exit("No experiments processed. Check config locations and input folders.")

    os.makedirs(results_dir, exist_ok=True)
    df_all = pd.DataFrame(all_rows).sort_values(["plot_label", "calendar_yr"])
    df_summary = pd.DataFrame(summary_rows).sort_values("Nae")

    xlsx_path = os.path.join(results_dir, "prr_all_experiments.xlsx")
    try:
        with pd.ExcelWriter(xlsx_path, engine="openpyxl") as writer:
            df_all.to_excel(writer, sheet_name="timeseries_all", index=False)
            df_summary.to_excel(writer, sheet_name="summary_by_experiment", index=False)
        print(f"\nSaved workbook: {xlsx_path}")
    except ImportError:
        all_csv = os.path.join(results_dir, "prr_all_experiments_timeseries.csv")
        summary_csv = os.path.join(results_dir, "prr_all_experiments_summary.csv")
        df_all.to_csv(all_csv, index=False)
        df_summary.to_csv(summary_csv, index=False)
        print(f"\nopenpyxl is not installed; saved CSVs: {all_csv}, {summary_csv}")

    plt.figure(figsize=(11, 6))
    for plot_label, grp in df_all.groupby("plot_label", sort=True):
        x = grp["model_yr"] if np.isfinite(grp["model_yr"]).any() else grp["calendar_yr"]
        plt.plot(x, grp["PRR"], marker="o", ms=2.5, lw=1, label=plot_label)
    plt.axhline(1.0, color="gray", ls="--", lw=0.8)
    plt.xlabel("Model year")
    plt.ylabel("PRR")
    plt.title("PRR through time by experiment")
    plt.grid(True, alpha=0.3)
    plt.legend(loc="best", fontsize=8, ncol=2)
    plt.tight_layout()
    time_plot_path = os.path.join(results_dir, "prr_vs_time_all_experiments.png")
    plt.savefig(time_plot_path, dpi=220, facecolor="white")
    plt.close()
    print(f"Saved plot: {time_plot_path}")

    _plot_summary(
        df_summary,
        "Nae",
        "Nae",
        os.path.join(results_dir, "prr_vs_nae_all_experiments.png"),
    )
    _plot_summary(
        df_summary,
        "slip_rate_mm_yr",
        "Slip rate (mm/yr)",
        os.path.join(results_dir, "prr_vs_slip_rate_all_experiments.png"),
    )


if __name__ == "__main__":
    main()
