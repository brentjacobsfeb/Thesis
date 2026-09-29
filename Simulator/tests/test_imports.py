"""Smoke tests for the simulator's initial component modules."""

import adc
import lna
import radar


def test_component_modules_import():
    assert adc.__name__ == "adc"
    assert lna.__name__ == "lna"
    assert radar.__name__ == "radar"