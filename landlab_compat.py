"""Landlab API differences across versions (constructor keyword names)."""
from __future__ import annotations

from inspect import signature
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from landlab import ModelGrid


def exponential_weatherer(grid: "ModelGrid", maximum_rate: float, decay_depth: float):
    """Build ``ExponentialWeatherer`` with kwargs that match the installed Landlab."""
    from landlab.components import ExponentialWeatherer

    sig = signature(ExponentialWeatherer.__init__)
    names = set(sig.parameters) - {"self", "cls"}

    kwargs = {}
    # Some Landlab versions use field-style names (double __ before max/decay).
    if "soil_production__maximum_rate" in names:
        kwargs["soil_production__maximum_rate"] = maximum_rate
    elif "soil_production_maximum_rate" in names:
        kwargs["soil_production_maximum_rate"] = maximum_rate
    elif "max_weathering_rate" in names:
        kwargs["max_weathering_rate"] = maximum_rate
    elif "maximum_weathering_rate" in names:
        kwargs["maximum_weathering_rate"] = maximum_rate
    else:
        raise TypeError(
            "ExponentialWeatherer has no recognized rate parameter; "
            f"got parameters: {sorted(names)}"
        )

    if "soil_production__decay_depth" in names:
        kwargs["soil_production__decay_depth"] = decay_depth
    elif "soil_production_decay_depth" in names:
        kwargs["soil_production_decay_depth"] = decay_depth
    else:
        raise TypeError(
            "ExponentialWeatherer has no recognized decay-depth parameter; "
            f"got parameters: {sorted(names)}"
        )

    return ExponentialWeatherer(grid, **kwargs)
