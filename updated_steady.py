#!/usr/bin/env python3
"""Build a clean steady-state topography for strike-slip experiments.

The runner is config-driven and hub-friendly. It saves:

- final steady-state grid as pickle and NetCDF
- copied final state in ``output/steady_state_files/``
- steady-state metrics CSV
- mean topography through time
- uplift/erosion and dz/dt diagnostics
- final topography
- final soil depth
- final slope-area plot
- optional paired soil-depth/topography evolution MP4, with temporary PNG frames deleted
"""
from __future__ import annotations

import argparse
import json
import pickle
from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace

import imageio.v2 as imageio
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.colors import LightSource
import numpy as np
import numpy_compat  # noqa: F401
import pandas as pd
import yaml
from landlab import RasterModelGrid, imshow_grid
from landlab.components import (
    DepthDependentTaylorDiffuser,
    PriorityFloodFlowRouter,
)
from landlab.components.space import SpaceLargeScaleEroder
from landlab.io.netcdf import write_netcdf

from landlab_compat import exponential_weatherer

ROOT = Path(__file__).resolve().parent
DEFAULT_CONFIG = ROOT / "config" / "steady_baseline_detachment.yaml"

NETCDF_FIELDS = [
    "bedrock__elevation",
    "drainage_area",
    "flood_status_code",
    "flow__link_to_receiver_node",
    "flow__receiver_node",
    "flow__receiver_proportions",
    "flow__upstream_node_order",
    "soil__depth",
    "soil_production__rate",
    "surface_water__discharge",
    "topographic__elevation",
    "topographic__steepest_slope",
    "water__unit_flux_in",
    "sediment__influx",
    "sediment__outflux",
    "sediment__flux",
]


def _section(data: dict, name: str) -> dict:
    section = data.get(name)
    if section is None:
        raise KeyError(f"Missing required config section: {name}")
    return section


def _geomorph_value(geomorphology: dict, new_key: str, old_key: str | None = None) -> float:
    if new_key in geomorphology:
        return float(geomorphology[new_key])
    if old_key is not None and old_key in geomorphology:
        return float(geomorphology[old_key])
    if old_key is None:
        raise KeyError(new_key)
    raise KeyError(f"{new_key} or {old_key}")


def load_config(path: Path) -> SimpleNamespace:
    with path.open("r") as f:
        data = yaml.safe_load(f)

    saving = _section(data, "saving")
    shape = _section(data, "shape")
    geomorphology = _section(data, "geomorphology")
    time = _section(data, "time")
    diagnostics = data.get("diagnostics", {})
    comments = data.get("comments", {})
    tectonics = data.get("tectonics", {})
    climate = data.get("climate", {})

    dxy = float(shape["dxy"])
    cfg = {
        "config_path": str(path),
        "raw_config": data,
        "model_name": str(saving["model_name"]),
        "alt_name": str(comments.get("alt_name", "")),
        "home_path": str(saving.get("home_path", ROOT)),
        "save_format": str(saving.get("output_filetype", "netcdf")),
        "frequency_output": float(saving.get("frequency_output", 10000.0)),
        "xmax": float(shape["xmax"]),
        "ymax": float(shape["ymax"]),
        "dxy": dxy,
        "ncols": int(round(float(shape["xmax"]) / dxy)),
        "nrows": int(round(float(shape["ymax"]) / dxy)),
        "H0": float(geomorphology["H0"]),
        "uplift_rate": float(geomorphology["uplift_rate"]),
        "Sc": float(geomorphology["Sc"]),
        "Hstar_d": _geomorph_value(geomorphology, "Hstar_d", "Hstar"),
        "Hstar_w": _geomorph_value(geomorphology, "Hstar_w", "Hstar"),
        "V0": float(geomorphology["V0"]),
        "taylor_nterms": int(geomorphology.get("taylor_nterms", 2)),
        "courant_factor": float(geomorphology.get("courant_factor", 0.1)),
        "P0": float(geomorphology["P0"]),
        "run_off": float(geomorphology["run_off"]),
        "K_sed": float(geomorphology["K_sed"]),
        "K_br": float(geomorphology["K_br"]),
        "F_f": float(geomorphology["F_f"]),
        "phi": float(geomorphology["phi"]),
        "Hstar_f": _geomorph_value(geomorphology, "Hstar_f", "H_star"),
        "Vs": float(geomorphology["Vs"]),
        "m_sp": float(geomorphology["m_sp"]),
        "n_sp": float(geomorphology["n_sp"]),
        "sp_crit_sed": float(geomorphology["sp_crit_sed"]),
        "sp_crit_br": float(geomorphology["sp_crit_br"]),
        "total_steady_time": float(time.get("total_steady_time", 1_000_000.0)),
        "dt_steady": float(time.get("dt_steady", 100.0)),
        "total_model_time": float(time.get("total_model_time", 0.0)),
        "dt_model": float(time.get("dt_model", 100.0)),
        "diagnostic_interval": float(diagnostics.get("diagnostic_interval", 5000.0)),
        "save_video": bool(diagnostics.get("save_video", True)),
        "video_interval": float(diagnostics.get("video_interval", 10000.0)),
        "video_fps": int(diagnostics.get("video_fps", 8)),
        "delete_video_frames": bool(diagnostics.get("delete_video_frames", True)),
        "random_seed": int(diagnostics.get("random_seed", 5000)),
        "initial_noise_m": float(diagnostics.get("initial_noise_m", 1.0)),
        "total_slip": float(tectonics.get("total_slip", 0.0)),
        "method": str(tectonics.get("method", "roll")),
        "slip_rate": float(tectonics.get("slip_rate", 0.0)),
        "fluvial_freq": float(climate.get("fluvial_freq", 100000.0)),
        "fluvial_len": float(climate.get("fluvial_len", 15000.0)),
    }
    cfg["D"] = cfg["Hstar_d"] * cfg["V0"]
    cfg["save_location"] = f"output/{cfg['model_name']}"
    return SimpleNamespace(**cfg)


def init_grid(config: SimpleNamespace) -> RasterModelGrid:
    mg = RasterModelGrid((config.nrows, config.ncols), config.dxy)
    mg.add_zeros("node", "topographic__elevation")

    rng = np.random.default_rng(config.random_seed)
    mg.at_node["topographic__elevation"] += (
        rng.random(mg.number_of_nodes) * config.initial_noise_m
    )

    mg.add_zeros("node", "soil__depth", clobber=True)
    mg.at_node["soil__depth"][mg.core_nodes] = config.H0

    mg.add_zeros("node", "bedrock__elevation", clobber=True)
    mg.at_node["bedrock__elevation"][:] = mg.at_node["topographic__elevation"]
    mg.at_node["topographic__elevation"][:] += mg.at_node["soil__depth"]
    mg.add_zeros("node", "soil_production__rate", clobber=True)

    mg.set_closed_boundaries_at_grid_edges(
        bottom_is_closed=False,
        left_is_closed=True,
        right_is_closed=True,
        top_is_closed=True,
    )
    return mg


def write_config(config: SimpleNamespace, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    serializable = {
        key: value
        for key, value in vars(config).items()
        if key != "raw_config"
    }
    with (output_dir / "config.json").open("w") as f:
        json.dump(serializable, f, indent=4, default=str)
    with (output_dir / "config_used.yaml").open("w") as f:
        yaml.safe_dump(config.raw_config, f, sort_keys=False)


def _shaded_image(
    field: np.ndarray,
    shade_source: np.ndarray,
    cmap_name: str,
    dxy: float,
    vmin: float | None = None,
    vmax: float | None = None,
) -> np.ndarray:
    if vmax is None:
        vmax = float(np.nanmax(field))
    if vmin is None:
        vmin = float(np.nanmin(field))
    if np.isclose(vmin, vmax):
        vmax = vmin + 1.0
    norm = colors.Normalize(
        vmin=vmin,
        vmax=vmax,
    )
    rgb = plt.get_cmap(cmap_name)(norm(field))
    hillshade = LightSource(azdeg=315, altdeg=45).hillshade(
        shade_source,
        vert_exag=1.0,
        dx=dxy,
        dy=dxy,
    )
    rgb[..., :3] *= 0.55 + 0.45 * hillshade[..., np.newaxis]
    return rgb


def make_soil_topography_frame(
    mg: RasterModelGrid,
    config: SimpleNamespace,
    frame_dir: Path,
    time: float,
) -> Path:
    frame_dir.mkdir(parents=True, exist_ok=True)
    soil = mg.at_node["soil__depth"].reshape((config.nrows, config.ncols))
    topo = mg.at_node["topographic__elevation"].reshape((config.nrows, config.ncols))
    extent = (0.0, config.xmax, 0.0, config.ymax)
    soil_vmin = float(np.nanmin(soil))
    soil_vmax = float(np.nanmax(soil))
    if np.isclose(soil_vmin, soil_vmax):
        soil_vmax = soil_vmin + 1.0
    soil_img = _shaded_image(soil, topo, "viridis", config.dxy, vmin=soil_vmin, vmax=soil_vmax)
    topo_img = _shaded_image(topo, topo, "terrain", config.dxy)

    fig, axes = plt.subplots(
        2,
        1,
        figsize=(8.5, 4.8),
        sharex=True,
        constrained_layout=True,
    )

    axes[0].imshow(
        soil_img,
        origin="lower",
        extent=extent,
        aspect="equal",
    )
    axes[0].set_ylabel("Y (m)")
    axes[0].set_title(f"{config.model_name}: soil depth, {time:,.0f} yr")
    soil_im = plt.cm.ScalarMappable(
        norm=colors.Normalize(vmin=soil_vmin, vmax=soil_vmax),
        cmap="viridis",
    )
    fig.colorbar(soil_im, ax=axes[0], label="Soil depth (m)")

    axes[1].imshow(
        topo_img,
        origin="lower",
        extent=extent,
        aspect="equal",
    )
    axes[1].set_xlabel("X (m)")
    axes[1].set_ylabel("Y (m)")
    axes[1].set_title("Topography")
    topo_im = plt.cm.ScalarMappable(
        norm=colors.Normalize(vmin=float(np.nanmin(topo)), vmax=float(np.nanmax(topo))),
        cmap="terrain",
    )
    fig.colorbar(topo_im, ax=axes[1], label="Elevation (m)")
    for ax in axes:
        ax.set_xlim(0.0, config.xmax)
        ax.set_ylim(0.0, config.ymax)

    frame_path = frame_dir / f"{config.model_name}_soil_topography_{int(round(time)):010d}.png"
    fig.savefig(frame_path, dpi=180, facecolor="white")
    plt.close(fig)
    return frame_path


def write_video(frame_paths: list[Path], video_path: Path, fps: int) -> Path | None:
    if not frame_paths:
        return None
    try:
        with imageio.get_writer(str(video_path), format="FFMPEG", fps=fps, macro_block_size=1) as writer:
            for frame_path in frame_paths:
                writer.append_data(imageio.imread(frame_path))
        print(f"Saved evolution video to {video_path}")
        return video_path
    except ImportError:
        gif_path = video_path.with_suffix(".gif")
        frames = [imageio.imread(frame_path) for frame_path in frame_paths]
        imageio.mimsave(gif_path, frames, duration=1 / max(fps, 1))
        print(
            "FFMPEG is not installed; saved GIF evolution instead. "
            f"Install imageio-ffmpeg for MP4 output. GIF: {gif_path}"
        )
        return gif_path


def save_final_grid(mg: RasterModelGrid, config: SimpleNamespace, output_dir: Path) -> None:
    steady_dir = ROOT / "output" / "steady_state_files"
    steady_dir.mkdir(parents=True, exist_ok=True)

    for path in [
        output_dir / f"final_state_{config.model_name}.pkl",
        steady_dir / f"final_state_{config.model_name}.pkl",
    ]:
        with path.open("wb") as f:
            pickle.dump(deepcopy(mg), f)
        print(f"Saved final steady-state pickle to {path}")

    mg_to_save = deepcopy(mg)
    for field_name in (
        "flow__receiver_node",
        "flow__link_to_receiver_node",
        "flow__upstream_node_order",
        "flood_status_code",
    ):
        if field_name in mg_to_save.at_node:
            mg_to_save.at_node[field_name] = mg_to_save.at_node[field_name].astype(np.int64)

    names = [name for name in NETCDF_FIELDS if name in mg_to_save.at_node]
    for path in [
        output_dir / f"final_state_{config.model_name}.nc",
        steady_dir / f"final_state_{config.model_name}.nc",
    ]:
        write_netcdf(path, mg_to_save, names=names, format="NETCDF3_64BIT")
        print(f"Saved final steady-state NetCDF to {path}")


def save_slope_area_plot(
    mg: RasterModelGrid,
    config: SimpleNamespace,
    output_dir: Path,
) -> None:
    area = np.asarray(mg.at_node["drainage_area"], dtype=float)
    slope = np.asarray(mg.at_node["topographic__steepest_slope"], dtype=float)
    core = mg.core_nodes
    mask = (
        np.isfinite(area[core])
        & np.isfinite(slope[core])
        & (area[core] > 0)
        & (slope[core] > 0)
    )
    area_core = area[core][mask]
    slope_core = slope[core][mask]

    fig, ax = plt.subplots(figsize=(7, 5), constrained_layout=True)
    ax.scatter(area_core, slope_core, s=3, c="black", alpha=0.35, linewidths=0)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Drainage area (m2)")
    ax.set_ylabel("Slope (m/m)")
    ax.set_title(f"{config.model_name}: final slope-area")
    ax.grid(True, which="both", alpha=0.25)
    fig.savefig(output_dir / "final_slope_area.png", dpi=220, facecolor="white")
    plt.close(fig)


def save_diagnostic_plots(
    mg: RasterModelGrid,
    metrics: pd.DataFrame,
    config: SimpleNamespace,
    output_dir: Path,
) -> None:
    fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
    ax.plot(metrics["time_yr"], metrics["mean_elevation_m"], lw=1.6)
    ax.set_xlabel("Time (yr)")
    ax.set_ylabel("Mean topographic elevation (m)")
    ax.set_title(f"{config.model_name}: mean topography through time")
    ax.grid(True, color="0.88", lw=0.7)
    fig.savefig(output_dir / "mean_topography_vs_time.png", dpi=220, facecolor="white")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
    ax.plot(metrics["time_yr"], metrics["uplift_rate_m_yr"], label="uplift rate", lw=1.5)
    ax.plot(
        metrics["time_yr"],
        metrics["erosion_equiv_rate_m_yr"],
        label="erosion-equivalent rate",
        lw=1.5,
    )
    ax.set_xlabel("Time (yr)")
    ax.set_ylabel("Rate (m/yr)")
    ax.set_title(f"{config.model_name}: uplift and erosion-equivalent rate")
    ax.grid(True, color="0.88", lw=0.7)
    ax.legend()
    fig.savefig(output_dir / "uplift_erosion_vs_time.png", dpi=220, facecolor="white")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
    ax.plot(metrics["time_yr"], metrics["mean_net_dz_rate_m_yr"], label="mean dz/dt", lw=1.5)
    ax.plot(metrics["time_yr"], metrics["mean_abs_dz_rate_m_yr"], label="mean abs dz/dt", lw=1.5)
    ax.axhline(0, color="0.25", lw=0.9, ls="--")
    ax.set_xlabel("Time (yr)")
    ax.set_ylabel("Rate (m/yr)")
    ax.set_title(f"{config.model_name}: topographic change rate")
    ax.grid(True, color="0.88", lw=0.7)
    ax.legend()
    fig.savefig(output_dir / "delta_z_vs_time.png", dpi=220, facecolor="white")
    plt.close(fig)

    topo = mg.at_node["topographic__elevation"].reshape((config.nrows, config.ncols))
    soil = mg.at_node["soil__depth"].reshape((config.nrows, config.ncols))
    extent = (0.0, config.xmax, 0.0, config.ymax)

    fig, ax = plt.subplots(figsize=(8, 3), constrained_layout=True)
    topo_img = _shaded_image(topo, topo, "terrain", config.dxy)
    ax.imshow(topo_img, origin="lower", extent=extent, aspect="equal")
    ax.set_xlabel("X (m)")
    ax.set_ylabel("Y (m)")
    ax.set_title(f"{config.model_name}: final steady-state topography")
    topo_im = plt.cm.ScalarMappable(
        norm=colors.Normalize(vmin=float(np.nanmin(topo)), vmax=float(np.nanmax(topo))),
        cmap="terrain",
    )
    fig.colorbar(topo_im, ax=ax, label="Elevation (m)")
    fig.savefig(output_dir / "final_topography.png", dpi=220, facecolor="white")
    plt.close(fig)

    soil_vmin = float(np.nanmin(soil))
    soil_vmax = float(np.nanmax(soil))
    if np.isclose(soil_vmin, soil_vmax):
        soil_vmax = soil_vmin + 1.0
    fig, ax = plt.subplots(figsize=(8, 3), constrained_layout=True)
    soil_img = _shaded_image(soil, topo, "viridis", config.dxy, vmin=soil_vmin, vmax=soil_vmax)
    ax.imshow(soil_img, origin="lower", extent=extent, aspect="equal")
    ax.set_xlabel("X (m)")
    ax.set_ylabel("Y (m)")
    ax.set_title(f"{config.model_name}: final soil depth")
    soil_im = plt.cm.ScalarMappable(
        norm=colors.Normalize(vmin=soil_vmin, vmax=soil_vmax),
        cmap="viridis",
    )
    fig.colorbar(soil_im, ax=ax, label="Soil depth (m)")
    fig.savefig(output_dir / "final_soil_depth.png", dpi=220, facecolor="white")
    plt.close(fig)

    save_slope_area_plot(mg, config, output_dir)


def run_steady(config: SimpleNamespace) -> RasterModelGrid:
    output_dir = ROOT / config.save_location
    frame_dir = output_dir / "video_frames"
    write_config(config, output_dir)

    mg = init_grid(config)
    rock = mg.at_node["bedrock__elevation"]
    z = mg.at_node["topographic__elevation"]
    soil = mg.at_node["soil__depth"]

    expweath = exponential_weatherer(mg, config.P0, config.Hstar_w)
    diffuser = DepthDependentTaylorDiffuser(
        mg,
        slope_crit=config.Sc,
        soil_transport_velocity=config.V0,
        soil_transport_decay_depth=config.Hstar_d,
        nterms=config.taylor_nterms,
        dynamic_dt=True,
        if_unstable="warn",
        courant_factor=config.courant_factor,
    )
    router = PriorityFloodFlowRouter(
        mg,
        flow_metric="D8",
        separate_hill_flow=False,
        hill_flow_metric="Quinn",
        update_hill_flow_instantaneous=False,
        suppress_out=True,
        runoff_rate=config.run_off,
    )
    eroder = SpaceLargeScaleEroder(
        mg,
        K_sed=config.K_sed,
        K_br=config.K_br,
        F_f=config.F_f,
        phi=config.phi,
        H_star=config.Hstar_f,
        v_s=config.Vs,
        m_sp=config.m_sp,
        n_sp=config.n_sp,
        sp_crit_sed=config.sp_crit_sed,
        sp_crit_br=config.sp_crit_br,
    )

    records = []
    frame_paths: list[Path] = []
    time = 0.0
    next_diagnostic = 0.0
    next_video = 0.0
    step = 0
    core = mg.core_nodes

    while time < config.total_steady_time:
        dt = min(config.dt_steady, config.total_steady_time - time)
        z_before = z.copy()

        rock[core] += config.uplift_rate * dt
        z[core] += config.uplift_rate * dt
        z_after_uplift = z.copy()

        expweath.calc_soil_prod_rate()
        diffuser.run_one_step(dt)
        router.run_one_step()
        eroder.run_one_step(dt)

        if not (
            np.all(np.isfinite(z))
            and np.all(np.isfinite(rock))
            and np.all(np.isfinite(soil))
        ):
            raise RuntimeError(
                "Non-finite values appeared in the steady-state grid. "
                "Try a smaller dt_steady, lower V0, lower initial_noise_m, or larger Sc."
            )

        time += dt
        step += 1

        if config.save_video and (time >= next_video or time >= config.total_steady_time):
            frame_paths.append(make_soil_topography_frame(mg, config, frame_dir, time))
            next_video += config.video_interval

        if time >= next_diagnostic or time >= config.total_steady_time:
            dz_rate = (z[core] - z_before[core]) / dt
            net_dz_rate = float(np.mean(dz_rate))
            mean_abs_dz_rate = float(np.mean(np.abs(dz_rate)))
            max_abs_dz_rate = float(np.max(np.abs(dz_rate)))
            process_lowering_rate = float(np.mean((z_after_uplift[core] - z[core]) / dt))
            erosion_equiv_rate = float(config.uplift_rate - net_dz_rate)
            mean_slope = float(np.mean(mg.at_node["topographic__steepest_slope"][core]))

            records.append(
                {
                    "time_yr": time,
                    "step": step,
                    "uplift_rate_m_yr": config.uplift_rate,
                    "mean_net_dz_rate_m_yr": net_dz_rate,
                    "mean_abs_dz_rate_m_yr": mean_abs_dz_rate,
                    "max_abs_dz_rate_m_yr": max_abs_dz_rate,
                    "process_lowering_rate_m_yr": process_lowering_rate,
                    "erosion_equiv_rate_m_yr": erosion_equiv_rate,
                    "mean_elevation_m": float(np.mean(z[core])),
                    "mean_soil_depth_m": float(np.mean(soil[core])),
                    "mean_drainage_area_m2": float(np.mean(mg.at_node["drainage_area"][core])),
                    "mean_slope": mean_slope,
                }
            )
            print(
                f"{time:,.0f} yr | mean topo={np.mean(z[core]):.3f} m | "
                f"erosion~{erosion_equiv_rate:.4g} m/yr | "
                f"mean dz/dt={net_dz_rate:.4g} m/yr | "
                f"mean abs dz/dt={mean_abs_dz_rate:.4g} m/yr"
            )
            next_diagnostic += config.diagnostic_interval

    metrics = pd.DataFrame(records)
    metrics.to_csv(output_dir / "steady_state_metrics.csv", index=False)

    save_final_grid(mg, config, output_dir)
    save_diagnostic_plots(mg, metrics, config, output_dir)

    if config.save_video:
        video_path = output_dir / f"{config.model_name}_soil_topography_evolution.mp4"
        write_video(frame_paths, video_path, fps=config.video_fps)
        if config.delete_video_frames:
            for frame_path in frame_paths:
                frame_path.unlink(missing_ok=True)
            try:
                frame_dir.rmdir()
            except OSError:
                pass
            print("Deleted temporary video PNG frames.")

    return mg


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--model-name", help="Override model_name from the config.")
    parser.add_argument("--total-steady-time", type=float, help="Override steady run duration.")
    parser.add_argument("--no-video", action="store_true", help="Do not save MP4/video frames.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    config = load_config(args.config)
    if args.model_name:
        config.model_name = args.model_name
        config.save_location = f"output/{config.model_name}"
    if args.total_steady_time is not None:
        config.total_steady_time = args.total_steady_time
    if args.no_video:
        config.save_video = False
    run_steady(config)


if __name__ == "__main__":
    main()
