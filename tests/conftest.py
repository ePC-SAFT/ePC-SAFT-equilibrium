from __future__ import annotations

import pytest


def pytest_addoption(parser: pytest.Parser) -> None:
    parser.addoption(
        "--held2-live",
        action="store_true",
        default=False,
        help="print observational HELD2.0 progress while the real controller runs",
    )


@pytest.fixture
def held2_live(request: pytest.FixtureRequest) -> bool:
    return bool(request.config.getoption("--held2-live"))
