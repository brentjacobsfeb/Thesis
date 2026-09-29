# Radar-Based Ground-Speed Sensor

This repository is the working area for a thesis project on validating the design of a radar-based ground-speed sensor. The Python simulator is in `Simulator/`; design literature and thesis planning are kept separately.

## Simulator status

The simulator is currently a scaffold. `radar.py`, `lna.py`, and `adc.py` are placeholders, and `main.py` only imports those modules. No radar, signal-chain, or ground-speed model has been implemented yet, so the current tests check imports rather than sensor behavior.

## Setup

## Layout

- `simulator/`: simulator entry point.
- `simulator/components/`: simulator components.
- `simulator/tests/`: simulator tests.
- `literature/`: local literature material; this directory remains excluded from Git by the existing ignore rule.