"""Patch NumPy 1.x aliases removed in NumPy 2.0 before importing Landlab.

Older Landlab wheels reference ``np.int`` etc. in C extensions. Import this module
once after ``import numpy`` and before ``import landlab``:

    import numpy as np
    import numpy_compat  # noqa: F401
"""
import numpy as np

if not hasattr(np, "int"):
    np.int = np.int64  # noqa: A001
if not hasattr(np, "float"):
    np.float = np.float64  # noqa: A001
if not hasattr(np, "complex"):
    np.complex = np.complex128  # noqa: A001
