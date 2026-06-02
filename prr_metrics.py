"""Profile Relief Ratio (PRR) metrics and Nae for strike-slip landscape models.

PRR_DT:
    Near-fault relief divided by farther-upstream relief. Profiles are placed
    at 10% and 50% of the fault-to-divide distance. Relief can be measured as
    full-profile relief or as maximum local relief in a moving window.

Profile placement
-----------------
Near-fault (10%) and far-fault (50%) of the distance from fault_row to
divide_row, measured in grid rows on the hanging-wall side.

Reference: Duvall & Tucker, JGR Earth Surface 2015.
"""
from __future__ import annotations

import numpy as np


def profile_row_indices(fault_row: int, divide_row: int) -> tuple[int, int]:
    """Row indices for 10% and 50% profiles between fault and divide.

    Returns (row_near, row_far).
    """
    span = divide_row - fault_row
    if span <= 0:
        raise ValueError(
            f"divide_row ({divide_row}) must be > fault_row ({fault_row})"
        )
    row_near = fault_row + round(0.10 * span)
    row_far = fault_row + round(0.50 * span)
    return int(row_near), int(row_far)


def compute_prr_dt(
    z_2d: np.ndarray,
    row_near: int,
    row_far: int,
    x_slice: slice | None = None,
) -> dict:
    """Duvall-Tucker PRR using full along-strike relief.

    Parameters
    ----------
    z_2d : (nrows, ncols) elevation array
    row_near, row_far : row indices for 10% and 50% profiles
    x_slice : column slice (default ``slice(1, -1)`` = interior columns)

    Returns
    -------
    dict with keys ``PRR_DT``, ``R_near``, ``R_far``.
    """
    if x_slice is None:
        x_slice = slice(1, -1)

    z_near = z_2d[row_near, x_slice].astype(float)
    z_far = z_2d[row_far, x_slice].astype(float)

    r_near = float(np.nanmax(z_near) - np.nanmin(z_near))
    r_far = float(np.nanmax(z_far) - np.nanmin(z_far))

    prr_dt = r_near / r_far if r_far > 0 else np.nan
    return {"PRR_DT": prr_dt, "R_near": r_near, "R_far": r_far}


def _max_local_relief(profile: np.ndarray, window_nodes: int) -> float:
    """Maximum local relief in a moving window along one profile."""
    profile = np.asarray(profile, dtype=float)
    if window_nodes < 2:
        raise ValueError("window_nodes must be at least 2")
    if window_nodes > profile.size:
        window_nodes = profile.size

    windows = np.lib.stride_tricks.sliding_window_view(profile, window_nodes)
    local_relief = np.nanmax(windows, axis=1) - np.nanmin(windows, axis=1)
    return float(np.nanmax(local_relief))


def compute_prr_dt_local(
    z_2d: np.ndarray,
    row_near: int,
    row_far: int,
    window_nodes: int,
    x_slice: slice | None = None,
) -> dict:
    """Duvall-Tucker PRR using maximum local relief in a moving window.

    ``PRR_DT_local = R_near_local / R_far_local``.
    """
    if x_slice is None:
        x_slice = slice(1, -1)

    z_near = z_2d[row_near, x_slice].astype(float)
    z_far = z_2d[row_far, x_slice].astype(float)

    r_near = _max_local_relief(z_near, window_nodes)
    r_far = _max_local_relief(z_far, window_nodes)
    prr_dt = r_near / r_far if r_far > 0 else np.nan

    return {
        "PRR_DT_local": prr_dt,
        "R_near_local": r_near,
        "R_far_local": r_far,
        "window_nodes": int(window_nodes),
    }


def compute_prr_tar(
    z_2d: np.ndarray,
    row_near: int,
    row_far: int,
    x_slice: slice | None = None,
    z_min_threshold: float = 1e-6,
) -> dict:
    """Per-node elevation ratio z_far/z_near, then mean/std/range.

    Parameters
    ----------
    z_2d : (nrows, ncols) elevation array
    row_near, row_far : row indices
    x_slice : column slice (default interior)
    z_min_threshold : ignore columns where abs(z_near) < threshold

    Returns
    -------
    dict with ``PRR_TAR_mean``, ``PRR_TAR_std``, ``PRR_TAR_range``,
    ``PRR_TAR_median``, ``n_valid``.
    """
    if x_slice is None:
        x_slice = slice(1, -1)

    z_near = z_2d[row_near, x_slice].astype(float)
    z_far = z_2d[row_far, x_slice].astype(float)

    valid = np.abs(z_near) >= z_min_threshold
    if not np.any(valid):
        return {
            "PRR_TAR_mean": np.nan,
            "PRR_TAR_std": np.nan,
            "PRR_TAR_range": np.nan,
            "PRR_TAR_median": np.nan,
            "n_valid": 0,
        }

    ratios = z_far[valid] / z_near[valid]
    return {
        "PRR_TAR_mean": float(np.nanmean(ratios)),
        "PRR_TAR_std": float(np.nanstd(ratios)),
        "PRR_TAR_range": float(np.nanmax(ratios) - np.nanmin(ratios)),
        "PRR_TAR_median": float(np.nanmedian(ratios)),
        "n_valid": int(np.sum(valid)),
    }


def compute_nae(slip_rate_mm_yr: float, K_br: float, D: float) -> float:
    """Nondimensional advection-erosion number.

    Nae = v^2 / (K_br * D)
    where v is slip rate in m/yr.
    """
    v = slip_rate_mm_yr / 1000.0  # mm/yr -> m/yr
    denominator = K_br * D
    if denominator == 0:
        return np.nan
    return (v ** 2) / denominator
