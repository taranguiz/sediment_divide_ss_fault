"""Profile Relief Ratio (PRR) metrics and Nae for strike-slip landscape models.


PRR_swath_profile:
    A near-fault swath is averaged into one
    strike-parallel profile, a farther-upstream swath is averaged into a second
    profile, and PRR is the ratio of their full-profile reliefs.

Profile placement
-----------------
Near-fault swath spans just upstream of the fault to 10% of the distance from
fault_row to divide_row. Far-upstream swath spans 40% to 50% of that distance.

Reference: Duvall & Tucker, JGR Earth Surface 2015.
Calculation based on the scripts built for Duvall, A. R., and G. E. Tucker (2015), 
Dynamic Ridges and Valleys in a Strike-Slip Environment, 
J. Geophys. Res. Earth Surf., 120, 2016–2026, doi:10.1002/ 2015JF003618.
"""
from __future__ import annotations

import numpy as np


def far_fault_swath_row_indices(
    fault_row: int,
    divide_row: int,
    lower_fraction: float = 0.40,
    upper_fraction: float = 0.50,
) -> tuple[int, int]:
    """Inclusive row bounds for an upstream swath between two fractions.

    Fractions are measured from ``fault_row`` to ``divide_row``. The returned
    pair is ``(row_start, row_end)`` and includes both endpoints.
    """
    span = divide_row - fault_row
    if span <= 0:
        raise ValueError(
            f"divide_row ({divide_row}) must be > fault_row ({fault_row})"
        )
    if not 0 <= lower_fraction <= upper_fraction <= 1:
        raise ValueError("swath fractions must satisfy 0 <= lower <= upper <= 1")

    row_start = fault_row + round(lower_fraction * span)
    row_end = fault_row + round(upper_fraction * span)
    return int(row_start), int(row_end)


def near_fault_swath_row_indices(
    fault_row: int,
    divide_row: int,
    upper_fraction: float = 0.10,
) -> tuple[int, int]:
    """Inclusive row bounds for the near-fault A swath.

    The default A swath spans from just upstream of the fault to 10% of the
    fault-to-divide distance.
    """
    span = divide_row - fault_row
    if span <= 0:
        raise ValueError(
            f"divide_row ({divide_row}) must be > fault_row ({fault_row})"
        )
    if not 0 <= upper_fraction <= 1:
        raise ValueError("upper_fraction must satisfy 0 <= upper_fraction <= 1")

    row_start = min(fault_row + 1, divide_row)
    row_end = fault_row + round(upper_fraction * span)
    row_end = max(row_start, min(row_end, divide_row))
    return int(row_start), int(row_end)


def x_slice_from_segment(
    ncols: int,
    dxy: float,
    segment_length_m: float | None = None,
    segment_start_m: float | None = None,
    edge_buffer_m: float = 0.0,
) -> slice:
    """Build a column slice for full-domain or fixed-length PRR segments.

    Parameters are in model distance units. If ``segment_length_m`` is omitted,
    the slice covers the interior domain after applying ``edge_buffer_m`` to
    both x edges. If ``segment_start_m`` is omitted for a fixed-length segment,
    the segment is centered in the buffered domain.
    """
    if ncols < 2:
        raise ValueError("ncols must be at least 2")
    if dxy <= 0:
        raise ValueError("dxy must be > 0")
    if edge_buffer_m < 0:
        raise ValueError("edge_buffer_m must be >= 0")

    first = int(np.ceil(edge_buffer_m / dxy))
    last_exclusive = ncols - int(np.ceil(edge_buffer_m / dxy))
    first = max(0, min(first, ncols - 1))
    last_exclusive = max(first + 1, min(last_exclusive, ncols))

    if segment_length_m is None:
        return slice(first, last_exclusive)
    if segment_length_m <= 0:
        raise ValueError("segment_length_m must be > 0")

    segment_nodes = max(2, int(round(segment_length_m / dxy)))
    available_nodes = last_exclusive - first
    if segment_nodes >= available_nodes:
        return slice(first, last_exclusive)

    if segment_start_m is None:
        start = first + (available_nodes - segment_nodes) // 2
    else:
        start = int(round(segment_start_m / dxy))
        start = max(first, min(start, last_exclusive - segment_nodes))
    return slice(start, start + segment_nodes)


def compute_prr_swath_profile(
    z_2d: np.ndarray,
    row_near_start: int,
    row_near_end: int,
    row_far_start: int,
    row_far_end: int,
    x_slice: slice | None = None,
) -> dict:
    """PRR from swath-averaged topographic profiles.

    Each swath is averaged across rows at every x column. Relief is then the
    full range of that averaged profile, matching the profile files read by
    ``dimensional_extractor.m`` from Duvall and Tucker, 2015.
    """
    if x_slice is None:
        x_slice = slice(1, -1)
    if row_near_end < row_near_start:
        raise ValueError("row_near_end must be >= row_near_start")
    if row_far_end < row_far_start:
        raise ValueError("row_far_end must be >= row_far_start")

    z_near = z_2d[row_near_start:row_near_end + 1, x_slice].astype(float)
    z_far = z_2d[row_far_start:row_far_end + 1, x_slice].astype(float)
    near_profile = np.nanmean(z_near, axis=0)
    far_profile = np.nanmean(z_far, axis=0)

    r_near = float(np.nanmax(near_profile) - np.nanmin(near_profile))
    r_far = float(np.nanmax(far_profile) - np.nanmin(far_profile))
    prr = r_near / r_far if r_far > 0 else np.nan

    return {
        "PRR_swath_profile": prr,
        "R_near_swath_profile": r_near,
        "R_far_swath_profile": r_far,
        "row_near_swath_start": int(row_near_start),
        "row_near_swath_end": int(row_near_end),
        "row_far_swath_start": int(row_far_start),
        "row_far_swath_end": int(row_far_end),
        "x_col_start": int(x_slice.start or 0),
        "x_col_end_exclusive": (
            int(x_slice.stop) if x_slice.stop is not None else int(z_2d.shape[1])
        ),
        "n_x_cols": int(near_profile.size),
        "near_profile": near_profile,
        "far_profile": far_profile,
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
