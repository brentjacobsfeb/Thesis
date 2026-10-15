# Radar-Based Ground-Speed Sensor

This repository is the working area for a thesis project on validating the design of a radar-based ground speed sensor.
The primary language is MATLAB.

## Simulator

The simulator models the functionality and non-idealities of the building blocks of the measurement system separately. The `main` function combines all the building blocks to a complete measurement chain for full system analysis.

A `Waveform` object carries the signal along the chain. Every component takes a waveform and returns a modified copy, so each stage can be kept and plotted:

```matlab
tx       = Waveform(samples, fs, fc, Label="Transmitted CW");
radiated = ant.transmit(tx);                    % radiated with the antenna pattern
echo     = ground.scatter(radiated, height);    % road echo back at the radar
rx       = ant.receive(echo);                   % signal at the antenna port
plot([tx, echo, rx])                            % time signal and spectrum of every stage
```

## Layout

- `simulator/`: simulator entry point.
- `simulator/components/`: These scipts model the functionality and non-idealities of the seperate components.
- `simulator/tests/`: These scripts conduct tests on every individual component.
- `simulator/components/main.m`: The complete measurement chain: transmits a CW carrier, scatters it off the moving road, receives the echo and estimates the ground speed from the Doppler centroid. Every set-up parameter is a name-value option. See `help main`.
- `simulator/tests/main_test.m`: Runs the complete measurement chain once with the default set-up, with all figures.
- `simulator/components/Waveform.m`: The signal passed between the components: complex baseband samples, sample rate and carrier frequency, with time-domain and spectrum plots. See `help Waveform`.
- `simulator/tests/WaveformTest.m`: Unit tests of the waveform.
- `simulator/components/RadarAntenna.m`: Antenna model: interpolated gain table, operating frequency and mounting on the vehicle; `transmit` and `receive` a waveform. See `help RadarAntenna`.
- `simulator/tests/RadarAntennaTest.m`: Unit tests of the antenna model.
- `simulator/components/GroundPlane.m`: Road model: constant-gamma backscatter (`sigma0 = gamma * sin(grazing)`) from a ground moving along the vehicle x axis; `scatter` turns a radiated waveform into the Doppler-shifted echo arriving back at the antenna. See `help GroundPlane`.
- `simulator/tests/GroundPlaneTest.m`: Unit tests of the road model.
- `simulator/run_tests.m`: Runs the tests.

## Running the simulator

From the `simulator/` folder in MATLAB:

```matlab
addpath components
main()                                              % default set-up, with figures
result = main(VehicleSpeed=30, Yaw=0, Plot=false);  % another set-up, no figures
result.EstimatedSpeed
```

## Running the tests

From the `simulator/` folder in MATLAB:

```matlab
run_tests()                                         % all tests
run_tests("RadarAntennaTest")                       % one test file
run_tests("RadarAntennaTest/taperLowersSidelobes")  % one test
```

From a shell: `matlab -batch "cd simulator; r = run_tests(); assertSuccess(r)"`.
