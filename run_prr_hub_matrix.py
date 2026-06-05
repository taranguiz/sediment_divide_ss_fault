#!/usr/bin/env python3
"""Run PRR-after-quake experiment variants from a committed run matrix.

This script is intended for the OpenEarthscape/JupyterHub workflow:

1. Clone the GitHub repo.
2. Upload required steady-state files to ``output/steady_state_files/``.
3. Run this script for one model at a time or for the enabled matrix rows.

Outputs are sparse NetCDF snapshots plus PRR/event tables. Topography PNG
frames and MP4 videos are off by default to keep disk use small.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
from pathlib import Path
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parent
DEFAULT_MATRIX = ROOT / "config" / "diffusion4_run_matrix.csv"
DEFAULT_TEMPLATES = ROOT / "config" / "diffusion4_family_templates.json"


def _read_matrix(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def _read_templates(path: Path) -> dict[str, dict]:
    with path.open() as f:
        return json.load(f)


def _sparse_output_frequency(total_model_time: float, dt_model: float, max_outputs: int) -> float:
    """Choose an interval that gives at most max_outputs saved model times."""
    if max_outputs < 2:
        return float(total_model_time + dt_model)
    interval = total_model_time / float(max_outputs - 1)
    steps = max(1, math.ceil(interval / dt_model))
    return float(steps * dt_model)


def _model_time_for_total_slip(total_slip: float, slip_rate_mm_yr: float, dt_model: float) -> float:
    """Return model duration needed to reach at least total_slip of slip."""
    if slip_rate_mm_yr <= 0:
        raise ValueError("slip_rate_mm_yr must be positive")
    slip_rate_m_yr = slip_rate_mm_yr / 1000.0
    raw_time = total_slip / slip_rate_m_yr
    steps = max(1, math.ceil(raw_time / dt_model))
    return float(steps * dt_model)


def _template_to_nested_config(model_name: str, template: dict, slip_rate: float, frequency: float) -> dict:
    return {
        "saving": {
            "model_name": model_name,
            "home_path": str(ROOT),
            "output_filetype": "netcdf",
            "frequency_output": frequency,
        },
        "comments": {"alt_name": template["alt_name"]},
        "shape": {
            "ymax": template["ymax"],
            "xmax": template["xmax"],
            "dxy": template["dxy"],
        },
        "geomorphology": {
            key: template[key]
            for key in [
                "H0",
                "uplift_rate",
                "Sc",
                "Hstar",
                "V0",
                "P0",
                "run_off",
                "K_sed",
                "K_br",
                "F_f",
                "phi",
                "H_star",
                "Vs",
                "m_sp",
                "n_sp",
                "sp_crit_sed",
                "sp_crit_br",
            ]
        },
        "tectonics": {
            "total_slip": template["total_slip"],
            "method": template["method"],
            "slip_rate": slip_rate,
        },
        "time": {
            "total_model_time": template["total_model_time"],
            "total_steady_time": template["total_steady_time"],
            "dt_steady": template["dt_steady"],
            "dt_model": template["dt_model"],
        },
        "climate": {
            "fluvial_freq": template["fluvial_freq"],
            "fluvial_len": template["fluvial_len"],
        },
    }


def _build_config(
    row: dict[str, str],
    template: dict,
    *,
    save_topo_plots: bool,
    delete_video_frames: bool,
    video_frame_mode: str,
) -> SimpleNamespace:
    model_name = row["model_name"]
    slip_rate = float(row["slip_rate_mm_yr"])
    total_slip = float(row["total_slip"])
    max_outputs = int(row["max_netcdf_outputs"])
    dt_model = float(template["dt_model"])
    total_model_time = _model_time_for_total_slip(total_slip, slip_rate, dt_model)

    data = dict(template)
    data["model_name"] = model_name
    data["plot_label"] = row["plot_label"]
    data["family"] = row["family"]
    data["alt_name"] = f'{row["plot_label"]}: {template["alt_name"]}'
    data["home_path"] = str(ROOT)
    data["save_location"] = f"output/{model_name}"
    data["save_format"] = "netcdf"
    data["slip_rate"] = slip_rate
    data["total_slip"] = total_slip
    data["total_model_time"] = total_model_time
    data["frequency_output"] = _sparse_output_frequency(
        total_model_time,
        dt_model,
        max_outputs,
    )
    data["initial_state_model_name"] = row["initial_state_model_name"]
    data["steady_state_file"] = row["steady_state_file"]
    data["max_netcdf_outputs"] = max_outputs
    data["sample_prr_at_quakes"] = True
    data["save_topo_plots"] = save_topo_plots
    data["topo_plot_frequency"] = data["frequency_output"]
    data["delete_video_frames"] = delete_video_frames
    data["video_frame_mode"] = video_frame_mode

    data["nrows"] = int(float(data["ymax"]) / float(data["dxy"]))
    data["ncols"] = int(float(data["xmax"]) / float(data["dxy"]))
    data["D"] = float(data["Hstar"]) * float(data["V0"])
    data["config"] = _template_to_nested_config(
        model_name,
        data,
        slip_rate,
        data["frequency_output"],
    )
    data["config"]["time"]["total_model_time"] = total_model_time
    data["config"]["tectonics"]["total_slip"] = total_slip

    return SimpleNamespace(**data)


def _select_rows(rows: list[dict[str, str]], labels: list[str] | None, include_disabled: bool) -> list[dict[str, str]]:
    if labels:
        wanted = set(labels)
        selected = [
            row for row in rows
            if row["plot_label"] in wanted or row["model_name"] in wanted
        ]
        missing = wanted.difference({row["plot_label"] for row in selected}).difference(
            {row["model_name"] for row in selected}
        )
        if missing:
            raise SystemExit(f"Unknown model label(s): {', '.join(sorted(missing))}")
        return selected

    if include_disabled:
        return rows
    return [row for row in rows if row["enabled"].strip().lower() == "yes"]


def _check_steady_state(row: dict[str, str]) -> Path:
    path = ROOT / "output" / "steady_state_files" / row["steady_state_file"]
    if not path.exists():
        raise FileNotFoundError(
            f"Missing steady-state file for {row['plot_label']}: {path}\n"
            "Upload/copy this file to output/steady_state_files/ before running."
        )
    return path


def _write_config_snapshot(config: SimpleNamespace) -> None:
    output_dir = ROOT / config.save_location
    output_dir.mkdir(parents=True, exist_ok=True)
    config_path = output_dir / "config.json"
    with config_path.open("w") as f:
        json.dump(vars(config), f, indent=4, default=str)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", type=Path, default=DEFAULT_MATRIX)
    parser.add_argument("--templates", type=Path, default=DEFAULT_TEMPLATES)
    parser.add_argument(
        "--label",
        action="append",
        help="Run only this plot label or model name. Can be passed multiple times.",
    )
    parser.add_argument(
        "--include-disabled",
        action="store_true",
        help="Run disabled matrix rows too, such as optional DT_10.",
    )
    parser.add_argument(
        "--with-topo-plots",
        action="store_true",
        help="Save topography PNG frames. Off by default to reduce hub storage.",
    )
    parser.add_argument(
        "--with-video",
        action="store_true",
        help="Save an MP4 evolution video. Implies --with-topo-plots.",
    )
    parser.add_argument(
        "--video-frame-mode",
        choices=["quake", "sparse"],
        default="quake",
        help="For --with-video, choose frames after quakes or sparse output intervals.",
    )
    parser.add_argument(
        "--keep-video-frames",
        action="store_true",
        help="Keep the temporary PNG frames used to build the MP4.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print selected runs and check steady-state files without running models.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    rows = _select_rows(_read_matrix(args.matrix), args.label, args.include_disabled)
    templates = _read_templates(args.templates)

    if not rows:
        raise SystemExit("No runs selected.")

    print(f"Selected {len(rows)} run(s). Repository root: {ROOT}")
    for row in rows:
        if row["family"] not in templates:
            raise SystemExit(f"No template found for family {row['family']}")
        steady_path = _check_steady_state(row)
        config = _build_config(
            row,
            templates[row["family"]],
            save_topo_plots=args.with_topo_plots or args.with_video,
            delete_video_frames=args.with_video and not args.keep_video_frames,
            video_frame_mode=args.video_frame_mode if args.with_video else "sparse",
        )
        print(
            f"{row['plot_label']} -> {config.model_name}: "
            f"slip={config.slip_rate:g} mm/yr, "
            f"total_slip={config.total_slip:g} m, "
            f"time={config.total_model_time:g} yr, "
            f"NetCDF every {config.frequency_output:g} yr, "
            f"steady={steady_path.name}"
        )
        if args.dry_run:
            continue

        from geomorph_dynamics_loop_trying_something import run_geomorf_loop

        _write_config_snapshot(config)
        writer = None
        if args.with_video:
            import imageio

            output_dir = ROOT / config.save_location
            output_dir.mkdir(parents=True, exist_ok=True)
            video_path = output_dir / f"{config.model_name}_evolution.mp4"
            writer = imageio.get_writer(video_path, fps=20)
            print(f"Saving evolution video to {video_path}")
        try:
            run_geomorf_loop(
                config,
                writer=writer,
                interactive_plots=False,
                save_outputs=True,
                sample_prr_at_quakes=True,
            )
        finally:
            if writer is not None:
                writer.close()


if __name__ == "__main__":
    os.chdir(ROOT)
    main()
