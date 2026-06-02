#!/usr/bin/env python3
"""
Compute Duvall-Tucker PRR (relief-based) on NetCDF experiment outputs.

This script can scan multiple experiment folders (each with .nc snapshots),
compute PRR_DT through time for each experiment, summarize mean/std PRR_DT
through time, and make combined PRR-vs-Nae and PRR-vs-slip-rate plots across
experiments.

PRR_DT = R_10 / R_50, matching the Duvall-Tucker near-fault / farther-upstream
definition. Relief is measured as max(z) - min(z) along each full
strike-parallel profile.

Usage:
    python test_prr_on_nc.py
"""
# %% Setup
import os
import sys
import glob
import json
import re

import numpy as np
import numpy_compat  # noqa: F401

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
}
ADDITIONAL_PRR_TIMESERIES = [
    {
        "experiment": "Duvall_Tucker_5_prr_full",
        "model_name": "Duvall_Tucker_5_prr_full",
        "plot_label": "DT_5_new",
        "csv_path": "/Users/taranguiz/Research/sediment_divide_ss_fault/output/Duvall_Tucker_5_prr_full/tabular_outputs/Duvall_Tucker_5_prr_full/Duvall_Tucker_5_prr_full_prr_at_quakes.csv",
        "config_path": "/Users/taranguiz/Research/sediment_divide_ss_fault/output/Duvall_Tucker_5_prr_full/config.json",
        "measurement_source": "quake_events",
    },
]

# %% Imports
from landlab.io.netcdf import read_netcdf

from prr_metrics import (
    profile_row_indices,
    compute_prr_dt,
    compute_nae,
)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


def _prompt_path(prompt: str, default_path: str) -> str:
    raw = input(f"{prompt} [{default_path}]: ").strip()
    return raw or default_path


def _extract_time_from_filename(path: str) -> int:
    name = os.path.basename(path)
    match = re.search(r"(\d+)\.nc$", name)
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
    """Return real NetCDF snapshots, ignoring macOS AppleDouble sidecar files."""
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
        if _netcdf_files(root):
            candidates = [root]
        else:
            candidates = [
                d for d in glob.glob(os.path.join(root, "*"))
                if os.path.isdir(d)
            ]
        for candidate in candidates:
            resolved = os.path.realpath(candidate)
            if resolved in seen:
                continue
            if _netcdf_files(candidate):
                exp_dirs.append(candidate)
                seen.add(resolved)
    return sorted(exp_dirs)


def _annotate_points(ax, df: pd.DataFrame, x_col: str, y_col: str) -> None:
    for _, row in df.iterrows():
        ax.annotate(
            row["plot_label"],
            (row[x_col], row[y_col]),
            xytext=(5, 5),
            textcoords="offset points",
            fontsize=8,
        )


def _plot_summary_ratio(
    df_summary: pd.DataFrame,
    x_col: str,
    x_label: str,
    y_col: str,
    yerr_col: str,
    y_label: str,
    title: str,
    output_path: str,
) -> None:
    fig, ax = plt.subplots(figsize=(8, 6))
    styles = {
        False: {
            "label": "uniform erodibility",
            "marker": "o",
            "color": "#2f6f9f",
        },
        True: {
            "label": "sediment effect (K_sed != K_br)",
            "marker": "s",
            "color": "#b84a39",
        },
    }
    for has_sediment_effect, grp in df_summary.groupby("has_sediment_effect", sort=True):
        style = styles[bool(has_sediment_effect)]
        ax.errorbar(
            grp[x_col],
            grp[y_col],
            yerr=grp[yerr_col],
            fmt=style["marker"],
            color=style["color"],
            ecolor=style["color"],
            capsize=4,
            elinewidth=1.2,
            linestyle="none",
            label=style["label"],
        )
    _annotate_points(ax, df_summary, x_col, y_col)
    ax.axhline(1.0, color="gray", ls="--", lw=0.8)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    fig.savefig(output_path, dpi=220, facecolor="white")
    print(f"Saved plot: {output_path}")


def _append_additional_prr_timeseries(
    all_rows: list[dict],
    summary_rows: list[dict],
) -> None:
    for extra in ADDITIONAL_PRR_TIMESERIES:
        csv_path = extra["csv_path"]
        config_path = extra["config_path"]
        if not os.path.exists(csv_path):
            print(f"Skipping {extra['plot_label']}: PRR CSV not found at {csv_path}")
            continue
        if not os.path.exists(config_path):
            print(f"Skipping {extra['plot_label']}: config not found at {config_path}")
            continue

        with open(config_path, "r") as f:
            cfg = json.load(f)
        exp_df = pd.read_csv(csv_path)
        if exp_df.empty:
            print(f"Skipping {extra['plot_label']}: PRR CSV is empty.")
            continue

        slip_rate = _slip_rate_from_config(cfg)
        k_sed = float(cfg["K_sed"])
        k_br = float(cfg["K_br"])
        d_coef = float(cfg["D"])
        dxy = float(cfg["dxy"])
        has_sediment_effect = not np.isclose(k_sed, k_br)
        nae = compute_nae(slip_rate, k_br, d_coef)
        first = exp_df.iloc[0]

        for _, row in exp_df.iterrows():
            all_rows.append(
                {
                    "experiment": extra["experiment"],
                    "model_name": extra["model_name"],
                    "plot_label": extra["plot_label"],
                    "snapshot": int(row["event_number"]),
                    "calendar_yr": float(row["calendar_year"]),
                    "model_yr": float(row["model_year"]),
                    "Nae": nae,
                    "slip_rate_mm_yr": slip_rate,
                    "K_sed": k_sed,
                    "K_br": k_br,
                    "has_sediment_effect": has_sediment_effect,
                    "D": d_coef,
                    "dxy": dxy,
                    "fault_row": int(row["fault_row"]),
                    "divide_row": int(row["divide_row"]),
                    "row_near": int(row["row_near"]),
                    "row_far": int(row["row_far"]),
                    "fault_y": float(row["fault_y"]),
                    "divide_y": float(row["divide_y"]),
                    "y_near": float(row["y_near"]),
                    "y_far": float(row["y_far"]),
                    "PRR_DT": float(row["PRR_DT"]),
                    "R_near": float(row["R_near"]),
                    "R_far": float(row["R_far"]),
                    "measurement_source": extra["measurement_source"],
                }
            )

        summary_rows.append(
            {
                "experiment": extra["experiment"],
                "model_name": extra["model_name"],
                "plot_label": extra["plot_label"],
                "Nae": nae,
                "slip_rate_mm_yr": slip_rate,
                "K_sed": k_sed,
                "K_br": k_br,
                "has_sediment_effect": has_sediment_effect,
                "D": d_coef,
                "dxy": dxy,
                "fault_row": int(first["fault_row"]),
                "divide_row": int(first["divide_row"]),
                "row_near": int(first["row_near"]),
                "row_far": int(first["row_far"]),
                "fault_y": float(first["fault_y"]),
                "divide_y": float(first["divide_y"]),
                "y_near": float(first["y_near"]),
                "y_far": float(first["y_far"]),
                "n_snapshots": len(exp_df),
                "measurement_source": extra["measurement_source"],
                "PRR_DT_mean_time": float(exp_df["PRR_DT"].mean()),
                "PRR_DT_std_time": float(exp_df["PRR_DT"].std()),
            }
        )
        print(
            f"\nAdded {extra['plot_label']}: {len(exp_df)} quake-event PRR measurements | "
            f"PRR mean={exp_df['PRR_DT'].mean():.4f} | std={exp_df['PRR_DT'].std():.4f}"
        )


def main():
    input_root_raw = _prompt_path(
        "Enter experiment output directories (comma-separated)",
        DEFAULT_INPUT_ROOTS,
    )
    input_roots = _split_input_roots(input_root_raw)
    local_config_root = _prompt_path("Enter local config root directory", DEFAULT_LOCAL_CONFIG_ROOT)
    results_dir = _prompt_path("Enter directory to save PRR results", DEFAULT_RESULTS_DIR)

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
            print(f"Skipping {exp_name}: config.json not found in experiment or local config root.")
            continue

        slip_rate = _slip_rate_from_config(cfg)
        k_sed = float(cfg["K_sed"])
        k_br = float(cfg["K_br"])
        d_coef = float(cfg["D"])
        has_sediment_effect = not np.isclose(k_sed, k_br)
        total_steady_time = int(cfg.get("total_steady_time", 0))
        model_name = _experiment_name_from_config(cfg, exp_dir)
        if model_name in DEFAULT_EXCLUDE_EXPERIMENTS:
            print(f"Skipping {model_name}: excluded from default PRR summary.")
            continue
        plot_label = EXPERIMENT_LABELS.get(model_name, model_name)
        nc_files = _netcdf_files(exp_dir)

        first_grid = read_netcdf(nc_files[0])
        nrows = int(first_grid.number_of_node_rows)
        ncols = int(first_grid.number_of_node_columns)
        dxy = float(first_grid.dx)
        fault_row = int(nrows / 3)
        divide_row = nrows - 1
        row_near, row_far = profile_row_indices(fault_row, divide_row)
        fault_y = float(first_grid.node_y[fault_row * ncols])
        divide_y = float(first_grid.node_y[divide_row * ncols])
        y_near = float(first_grid.node_y[row_near * ncols])
        y_far = float(first_grid.node_y[row_far * ncols])
        nae = compute_nae(slip_rate, k_br, d_coef)

        print(
            f"\n{exp_name}: {len(nc_files)} snapshots | rows {row_near}/{row_far} | "
            f"y={y_near:.1f}/{y_far:.1f} m | grid={nrows}x{ncols} | "
            f"Nae={nae:.4f} | config={cfg['_config_path']}"
        )

        exp_rows = []
        for i, nc_path in enumerate(nc_files):
            calendar_yr = _extract_time_from_filename(nc_path)
            model_yr = (calendar_yr - total_steady_time) if calendar_yr >= 0 else np.nan

            mg = read_netcdf(nc_path)
            file_nrows = int(mg.number_of_node_rows)
            file_ncols = int(mg.number_of_node_columns)
            if (file_nrows, file_ncols) != (nrows, ncols):
                raise ValueError(
                    f"{nc_path} has grid {file_nrows}x{file_ncols}, expected {nrows}x{ncols}"
                )
            z_2d = mg.at_node["topographic__elevation"].reshape((file_nrows, file_ncols))
            dt = compute_prr_dt(z_2d, row_near, row_far)

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
                "dxy": dxy,
                "fault_row": fault_row,
                "divide_row": divide_row,
                "row_near": row_near,
                "row_far": row_far,
                "fault_y": fault_y,
                "divide_y": divide_y,
                "y_near": y_near,
                "y_far": y_far,
                "PRR_DT": dt["PRR_DT"],
                "R_near": dt["R_near"],
                "R_far": dt["R_far"],
                "measurement_source": "netcdf_snapshots",
            }
            exp_rows.append(row)
            all_rows.append(row)

        if not exp_rows:
            continue

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
                "row_near": row_near,
                "row_far": row_far,
                "fault_y": fault_y,
                "divide_y": divide_y,
                "y_near": y_near,
                "y_far": y_far,
                "n_snapshots": len(exp_df),
                "measurement_source": "netcdf_snapshots",
                "PRR_DT_mean_time": float(exp_df["PRR_DT"].mean()),
                "PRR_DT_std_time": float(exp_df["PRR_DT"].std()),
            }
        )

    _append_additional_prr_timeseries(all_rows, summary_rows)

    if not all_rows:
        sys.exit("No experiments processed. Check config locations and input folders.")

    os.makedirs(results_dir, exist_ok=True)
    df_all = pd.DataFrame(all_rows).sort_values(["plot_label", "calendar_yr"])
    df_summary = pd.DataFrame(summary_rows).sort_values("Nae")

    xlsx_path = os.path.join(results_dir, "prr_dt_all_experiments.xlsx")
    with pd.ExcelWriter(xlsx_path, engine="openpyxl") as writer:
        df_all.to_excel(writer, sheet_name="timeseries_all", index=False)
        df_summary.to_excel(writer, sheet_name="summary_by_experiment", index=False)
    print(f"\nSaved workbook: {xlsx_path}")

    # PRR through time for each experiment
    plt.figure(figsize=(11, 6))
    for plot_label, grp in df_all.groupby("plot_label", sort=True):
        x = grp["model_yr"] if np.isfinite(grp["model_yr"]).any() else grp["calendar_yr"]
        plt.plot(x, grp["PRR_DT"], marker="o", ms=2.5, lw=1, label=plot_label)
    plt.axhline(1.0, color="gray", ls="--", lw=0.8)
    plt.xlabel("Model year")
    plt.ylabel("PRR_DT (R_10 / R_50)")
    plt.title("PRR_DT through time by experiment")
    plt.grid(True, alpha=0.3)
    plt.legend(loc="best", fontsize=8, ncol=2)
    plt.tight_layout()
    time_plot_path = os.path.join(results_dir, "prr_dt_vs_time_all_experiments.png")
    plt.savefig(time_plot_path, dpi=220, facecolor="white")
    print(f"Saved plot: {time_plot_path}")

    _plot_summary_ratio(
        df_summary,
        "Nae",
        "Nae",
        "PRR_DT_mean_time",
        "PRR_DT_std_time",
        "PRR_DT mean through time (+/- 1 std)",
        "PRR_DT vs Nae across experiments",
        os.path.join(results_dir, "prr_dt_vs_nae_all_experiments.png"),
    )
    _plot_summary_ratio(
        df_summary,
        "slip_rate_mm_yr",
        "Slip rate (mm/yr)",
        "PRR_DT_mean_time",
        "PRR_DT_std_time",
        "PRR_DT mean through time (+/- 1 std)",
        "PRR_DT vs slip rate across experiments",
        os.path.join(results_dir, "prr_dt_vs_slip_rate_all_experiments.png"),
    )

    plt.close("all")


if __name__ == "__main__":
    main()
