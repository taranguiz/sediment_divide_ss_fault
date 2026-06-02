#!/usr/bin/env python3
"""Run a local Duvall_Tucker_5-style test with PRR sampled after earthquakes.

By default this is a short smoke test. Set ``RUN_FULL=1`` in the environment
to run the full 100 kyr Duvall_Tucker_5 duration.
"""
from __future__ import annotations

import json
import math
import os
from types import SimpleNamespace

import imageio

from geomorph_dynamics_loop_trying_something import run_geomorf_loop


ROOT = "/Users/taranguiz/Research/sediment_divide_ss_fault"
SOURCE_CONFIG = os.path.join(ROOT, "output", "Duvall_Tucker_5", "config.json")
SMOKE_OUTPUT_MODEL_NAME = "Duvall_Tucker_5_prr_test"
FULL_OUTPUT_MODEL_NAME = "Duvall_Tucker_5_prr_full"
MAX_NETCDF_OUTPUTS = 10


def _load_config() -> SimpleNamespace:
    with open(SOURCE_CONFIG, "r") as f:
        data = json.load(f)

    config = SimpleNamespace(**data)
    config.config = data.get("config", {})
    config.initial_state_model_name = "Duvall_Tucker_5"
    is_full_run = os.environ.get("RUN_FULL") == "1"
    config.model_name = FULL_OUTPUT_MODEL_NAME if is_full_run else SMOKE_OUTPUT_MODEL_NAME
    config.alt_name = (
        "Duvall_Tucker_5_PRR_after_quakes_full"
        if is_full_run
        else "Duvall_Tucker_5_PRR_after_quakes_test"
    )
    config.home_path = ROOT
    config.save_location = f"output/{config.model_name}"
    config.save_format = "netcdf"

    if not is_full_run:
        config.total_model_time = 2000.0

    config.frequency_output = _sparse_output_frequency(
        config.total_model_time,
        config.dt_model,
        MAX_NETCDF_OUTPUTS,
    )
    return config


def _sparse_output_frequency(total_model_time, dt_model, max_outputs):
    """Choose an output interval that yields at most max_outputs unique grids."""
    if max_outputs < 2:
        return float(total_model_time + dt_model)
    interval = total_model_time / float(max_outputs - 1)
    steps = max(1, math.ceil(interval / dt_model))
    return float(steps * dt_model)


def main():
    config = _load_config()
    output_dir = os.path.join(config.home_path, config.save_location)
    os.makedirs(output_dir, exist_ok=True)
    with open(os.path.join(output_dir, "config.json"), "w") as f:
        json.dump(vars(config), f, indent=4, default=str)

    print(
        f"Running {config.model_name}: total_model_time={config.total_model_time:g}, "
        f"frequency_output={config.frequency_output:g}, "
        "PRR sampled after every earthquake."
    )
    video_path = os.path.join(output_dir, f"{config.model_name}_evolution.mp4")
    writer = imageio.get_writer(video_path, fps=20)
    try:
        run_geomorf_loop(
            config,
            writer=writer,
            interactive_plots=False,
            save_outputs=True,
            sample_prr_at_quakes=True,
        )
    finally:
        writer.close()
    print(f"Saved evolution video to {video_path}")


if __name__ == "__main__":
    main()
