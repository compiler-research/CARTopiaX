"""Unit tests for the behavior system.  No BioDynaMo required."""

import ctypes

import pytest

from cartopiaX.behavior import (
    HypoxiaApoptosis,
    PyBehavior,
    _PyGILState_Ensure,
    _PyGILState_Release,
    _registry,
    get_behavior,
    register_behavior,
)


# ── Registry ──────────────────────────────────────────────────────────────────

def test_hypoxia_registered_at_import():
    assert "HypoxiaApoptosis" in _registry


def test_get_behavior_returns_correct_class():
    assert get_behavior("HypoxiaApoptosis") is HypoxiaApoptosis


def test_get_behavior_raises_for_unknown_name():
    with pytest.raises(KeyError, match="no_such"):
        get_behavior("no_such")


def test_register_custom_and_retrieve():
    @register_behavior
    class _Tmp(PyBehavior):
        def run(self, agent):
            pass

    try:
        assert get_behavior("_Tmp") is _Tmp
    finally:
        del _registry["_Tmp"]


def test_register_overwrites_existing():
    @register_behavior
    class _Ow(PyBehavior):
        def run(self, agent):
            pass

    original = _registry.get("_Ow")

    @register_behavior
    class _Ow(PyBehavior):  # noqa: F811
        def run(self, agent):
            pass

    try:
        assert _registry["_Ow"] is not original
    finally:
        del _registry["_Ow"]


# ── PyBehavior GIL handles ────────────────────────────────────────────────────

def test_gil_handles_are_resolved():
    # If ctypes.pythonapi could not resolve these names the module would have
    # failed to import, but we assert the types explicitly for clarity.
    assert callable(_PyGILState_Ensure)
    assert callable(_PyGILState_Release)


def test_run_dispatches_to_run_method():
    """PyBehavior.Run() must call self.run() exactly once."""

    class _Counter(PyBehavior):
        calls = 0

        def run(self, agent):
            self.calls += 1

    b = _Counter()
    b.Run(object())
    assert b.calls == 1


def test_run_releases_gil_on_exception():
    """GIL must be released even when run() raises."""

    class _Boom(PyBehavior):
        def run(self, agent):
            raise ValueError("boom")

    b = _Boom()
    with pytest.raises(ValueError, match="boom"):
        b.Run(object())
    # If we reach here the GIL was released (otherwise the interpreter would hang).


# ── HypoxiaApoptosis ──────────────────────────────────────────────────────────

def test_default_threshold():
    assert HypoxiaApoptosis().threshold == 5.0


def test_custom_threshold():
    assert HypoxiaApoptosis(threshold=2.5).threshold == 2.5


class _MockGrid:
    def __init__(self, oxygen):
        self._oxygen = oxygen

    def GetValue(self, pos):
        return self._oxygen


class _MockAgent:
    def __init__(self, oxygen, has_grid=True):
        self._grid = _MockGrid(oxygen) if has_grid else None
        self.apoptosis_called = False

    def GetOxygenDiffusionGrid(self):
        return self._grid

    def GetPosition(self):
        return (0.0, 0.0, 0.0)

    def StartApoptosis(self):
        self.apoptosis_called = True


def test_no_apoptosis_above_threshold():
    agent = _MockAgent(oxygen=20.0)
    HypoxiaApoptosis(threshold=5.0).run(agent)
    assert not agent.apoptosis_called


def test_no_apoptosis_at_threshold():
    agent = _MockAgent(oxygen=5.0)
    HypoxiaApoptosis(threshold=5.0).run(agent)
    assert not agent.apoptosis_called


def test_apoptosis_triggered_below_threshold():
    agent = _MockAgent(oxygen=1.0)
    HypoxiaApoptosis(threshold=5.0).run(agent)
    assert agent.apoptosis_called


def test_no_error_when_grid_is_none():
    agent = _MockAgent(oxygen=0.0, has_grid=False)
    HypoxiaApoptosis().run(agent)
    assert not agent.apoptosis_called


def test_run_via_dispatch_triggers_apoptosis():
    """Full Run() path (GIL acquire → run() → GIL release) must also work."""
    agent = _MockAgent(oxygen=0.5)
    HypoxiaApoptosis(threshold=5.0).Run(agent)
    assert agent.apoptosis_called
