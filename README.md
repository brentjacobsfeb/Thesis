# Radar-Based Ground-Speed Sensor

This repository is the working area for a thesis project on validating the design of a radar-based ground speed sensor.
The primary language is MATLAB.

## Simulator

The simulator models the functionality and non-idealities of the building blocks of the measurement system separately. The main script combines all the building blocks to a complete measurement chain for full system analysis.

## Layout

- `simulator/`: simulator entry point.
- `simulator/components/`: These scipts model the functionality and non-idealities of the seperate components.
- `simulator/tests/`: These scripts conduct tests on every individual component.
- `simulator/components/main.m`: A complete measurement chain function.
- `simulator/components/main_test.m`: A complete measurement chain tester.
- `simulator/components/RadarAntenna.m`: Antenna model: interpolated gain table, operating frequency and mounting on the vehicle. See `help RadarAntenna`.
- `simulator/tests/RadarAntennaTest.m`: Unit tests of the antenna model.
- `simulator/components/GroundPlane.m`: Road model: constant-gamma backscatter (`sigma0 = gamma * sin(grazing)`) from a ground moving along the vehicle x axis, producing the Doppler-shifted echo the antenna receives. See `help GroundPlane`.
- `simulator/tests/GroundPlaneTest.m`: Unit tests of the road model.
- `simulator/demo.m`: Examples of the antenna model (table from a formula, ground footprint, ground echo and Doppler speed estimate, importing an Antenna Toolbox pattern).
- `simulator/run_tests.m`: Runs the tests.

## Running the tests

From the `simulator/` folder in MATLAB:

```matlab
run_tests()                                         % all tests
run_tests("RadarAntennaTest")                       % one test file
run_tests("RadarAntennaTest/taperLowersSidelobes")  % one test
```

From a shell: `matlab -batch "cd simulator; r = run_tests(); assertSuccess(r)"`.
