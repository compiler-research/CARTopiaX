"""cppyy bootstrap: load BioDynaMo and CARTopiaX shared libraries.

Lines to update for your installation
--------------------------------------
BDM_INSTALL_DIR   – root of the BioDynaMo install tree (contains lib/ and include/).
                    Override via the BDM_INSTALL_DIR environment variable.

CARTOPIX_LIB_DIR  – directory that contains libCARTopiaX.so.
                    BioDynaMo's bdm_add_executable macro builds a shared library
                    alongside the executable; this is typically <build>/lib/ or
                    the build root itself.
                    Override via the CARTOPIX_LIB_DIR environment variable.

CARTOPIX_SRC_DIR  – directory that contains the CARTopiaX .h headers
                    (hyperparams.h, tumor_cell.h, cart_cell.h, cart_tumor.h).
                    Override via the CARTOPIX_SRC_DIR environment variable.
"""

from __future__ import annotations

import os

# ── UPDATE THESE THREE PATHS FOR YOUR INSTALLATION ──────────────────────────
BDM_INSTALL_DIR  = os.environ.get("BDM_INSTALL_DIR",  "/opt/biodynamo")
CARTOPIX_LIB_DIR = os.environ.get("CARTOPIX_LIB_DIR", "/opt/cartopiaX/build")
CARTOPIX_SRC_DIR = os.environ.get("CARTOPIX_SRC_DIR", "/opt/cartopiaX/src")
# ────────────────────────────────────────────────────────────────────────────

_loaded = False


def load() -> None:
    """Load both shared libraries and parse headers into Cling.  Safe to call repeatedly."""
    global _loaded
    if _loaded:
        return

    import cppyy  # noqa: PLC0415 – intentionally deferred

    bdm_lib_dir = os.path.join(BDM_INSTALL_DIR, "lib")
    bdm_inc_dir = os.path.join(BDM_INSTALL_DIR, "include")

    cppyy.add_library_path(bdm_lib_dir)
    cppyy.add_library_path(CARTOPIX_LIB_DIR)

    # BioDynaMo must be loaded before the project library because CARTopiaX
    # links against it; the dynamic linker resolves symbols in load order.
    cppyy.load_library("libbdm")
    cppyy.load_library("libCARTopiaX")

    cppyy.add_include_path(bdm_inc_dir)
    # BioDynaMo installs its own headers under include/bdm/ in some builds.
    cppyy.add_include_path(os.path.join(bdm_inc_dir, "bdm"))
    cppyy.add_include_path(CARTOPIX_SRC_DIR)

    # biodynamo.h is the top-level convenience header; it pulls in simulation,
    # scheduler, resource manager, etc.  Project headers follow.
    cppyy.include("biodynamo.h")
    cppyy.include("hyperparams.h")
    cppyy.include("tumor_cell.h")
    cppyy.include("cart_cell.h")
    cppyy.include("cart_tumor.h")

    _loaded = True


def require() -> None:
    """Call load(); re-raise as ImportError with actionable guidance on failure."""
    try:
        load()
    except Exception as exc:
        raise ImportError(
            f"Could not load BioDynaMo/CARTopiaX libraries: {exc}\n"
            "Set these environment variables to your install paths:\n"
            "  BDM_INSTALL_DIR   – BioDynaMo install root\n"
            "  CARTOPIX_LIB_DIR  – directory containing libCARTopiaX.so\n"
            "  CARTOPIX_SRC_DIR  – directory containing the CARTopiaX headers"
        ) from exc


def cpp_treatment_map(py_dict: dict[int, int]):
    """Convert a Python {int: int} dict to a bdm::SimParam treatment std::map.

    Ownership: the returned map is a cppyy-managed C++ value; Python owns it
    until it is copied into a SimParam field.
    """
    require()
    import cppyy  # noqa: PLC0415

    m = cppyy.gbl.std.map[int, int]()
    for day, dose in py_dict.items():
        m[day] = dose
    return m
