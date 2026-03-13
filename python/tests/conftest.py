import pytest


def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "integration: marks tests that require BioDynaMo and CARTopiaX "
        "shared libraries to be installed and configured via environment variables.",
    )
