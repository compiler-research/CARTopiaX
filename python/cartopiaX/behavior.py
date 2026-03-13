"""Behavior base class, registration system, and example concrete behavior.

GIL handling
------------
BioDynaMo dispatches Behavior::Run() from OpenMP worker threads.  CPython's
GIL is not held on entry to any thread other than the main thread, so it must
be explicitly acquired before touching any Python object.

The acquire/release is inline in PyBehavior.Run() using raw ctypes handles to
PyGILState_Ensure / PyGILState_Release.  It is NOT hidden inside a decorator
so that a reviewer can audit the acquire/release pair directly.

Ownership
---------
PyBehavior instances passed to a C++ agent via AddBehavior() are owned by
BioDynaMo after that call.  Do not retain a Python reference and call methods
on the object afterwards — the C++ side may free it at any point.

Instances that have NOT been passed to AddBehavior() are owned by Python as
usual (e.g., in unit tests calling run() directly).
"""

from __future__ import annotations

import ctypes
from abc import ABC, abstractmethod
from typing import Any, Dict, Type

# ── Behavior registry ─────────────────────────────────────────────────────────

_registry: Dict[str, Type["PyBehavior"]] = {}


def register_behavior(cls: Type["PyBehavior"]) -> Type["PyBehavior"]:
    """Decorator: store a PyBehavior subclass in the module-level registry."""
    _registry[cls.__name__] = cls
    return cls


def get_behavior(name: str) -> Type["PyBehavior"]:
    """Look up a registered behavior class by name."""
    if name not in _registry:
        raise KeyError(
            f"No behavior registered under '{name}'; "
            f"registered names: {sorted(_registry)}"
        )
    return _registry[name]


# ── Raw GIL handles ───────────────────────────────────────────────────────────
# Resolved once at module import so Run() dispatch has no per-call overhead
# from attribute lookup on ctypes.pythonapi.

_PyGILState_Ensure = ctypes.pythonapi.PyGILState_Ensure
_PyGILState_Release = ctypes.pythonapi.PyGILState_Release
_PyGILState_Ensure.restype = ctypes.c_int
_PyGILState_Release.argtypes = [ctypes.c_int]


# ── Base class ────────────────────────────────────────────────────────────────

class PyBehavior(ABC):
    """Python base for CARTopiaX behaviors.

    Subclass and override run().  The C++-facing Run() acquires the GIL before
    delegating to run() and releases it after — unconditionally, even on
    exception.

    When BioDynaMo is available, use make_cpp_behavior_base() to obtain a
    version of this class that also inherits from cppyy's bdm::Behavior, which
    is required for AddBehavior() to accept the instance.
    """

    def Run(self, agent: Any) -> None:
        # Called from a BioDynaMo OpenMP thread.  The GIL is not held on entry.
        gil_state = _PyGILState_Ensure()
        try:
            self.run(agent)
        finally:
            # Release even if run() raises so we never leave the GIL leaked.
            _PyGILState_Release(gil_state)

    @abstractmethod
    def run(self, agent: Any) -> None:
        """Override in subclasses.  `agent` is a cppyy bdm::Agent* (or mock)."""


# ── cppyy-backed variant ──────────────────────────────────────────────────────

def make_cpp_behavior_base():
    """Return a class inheriting from both PyBehavior and cppyy's bdm::Behavior.

    Call this lazily (only when BioDynaMo is available).  The returned class
    can be used as a base so instances are accepted by bdm::Agent::AddBehavior.

    The Run() override is repeated explicitly here rather than relying on MRO
    so the GIL acquire/release pair remains unconditionally visible in this
    class without digging through parent classes.
    """
    from . import _bootstrap  # noqa: PLC0415
    import cppyy  # noqa: PLC0415

    _bootstrap.require()

    class _CppPyBehavior(cppyy.gbl.bdm.Behavior, PyBehavior):
        """cppyy-backed PyBehavior; use this base when adding behaviors to C++ agents.

        Structural scaffolding — the bdm::Behavior base class and the MRO it
        implies are untested without a real BioDynaMo install.  No current test
        exercises this class; it is not collected by the unit-test suite and the
        integration tests do not call make_cpp_behavior_base() either.  The
        inheritance is present so the plumbing is in place for a future phase
        that wires Python behaviors into BioDynaMo's per-step scheduler.
        """

        def Run(self, agent):
            # Identical to PyBehavior.Run() — duplicated so the GIL pair is
            # plainly visible here without needing to trace MRO.
            gil_state = _PyGILState_Ensure()
            try:
                self.run(agent)
            finally:
                _PyGILState_Release(gil_state)

        @abstractmethod
        def run(self, agent):
            pass

    return _CppPyBehavior


# ── Concrete example behavior ─────────────────────────────────────────────────

@register_behavior
class HypoxiaApoptosis(PyBehavior):
    """Trigger apoptosis on a TumorCell when local oxygen falls below a threshold.

    This demonstrates the pattern for Python-defined behaviors: read a cell
    property via the cppyy-wrapped getter, make a decision, call a mutating
    method if warranted.
    """

    def __init__(self, threshold: float = 5.0) -> None:
        """
        Args:
            threshold: Oxygen level in mmHg below which apoptosis is triggered.
                       Matches oxygen_limit_for_necrosis (5.0 mmHg) by default.
        """
        self.threshold = threshold

    def run(self, agent: Any) -> None:
        # GetOxygenDiffusionGrid() is defined on both TumorCell and CarTCell.
        # We do not attempt a dynamic cast here; if the agent lacks the method
        # (e.g., a mock without a grid), we bail silently.
        grid = agent.GetOxygenDiffusionGrid()
        if grid is None:
            return
        oxygen = grid.GetValue(agent.GetPosition())
        if oxygen < self.threshold:
            # StartApoptosis() is defined on TumorCell (cart_cell.h does not
            # expose it); cppyy will raise AttributeError at runtime if this
            # is called on a CarTCell, which is the correct failure mode.
            agent.StartApoptosis()
