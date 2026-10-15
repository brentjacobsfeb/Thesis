close all;
clear;
addpath(fullfile(fileparts(mfilename("fullpath")), "components"));

%% define antenna parameters
fc = 24.125e9;          % Carrier frequency [Hz]
mountHeight = 0.5;      % Antenna height above the road [m]
depression = 45;        % Boresight angle below the horizon [deg]
yaw = 180;              %
roll = 0;               %
az = -180:1:180;
el = (-90:1:90).';
gainDb = 20 - 10 * log10(2) * 4 * ((az / 12).^2 + (el / 20).^2);    % Gain table [dB]
phaseDeg = zeros(size(gainDb));                                     % Phase table [deg]

%% create antenna
ant = RadarAntenna(az, el, gainDb, phaseDeg, fc, depression, yaw, roll);
ant.plotCuts(-90:0.1:90, 60);
ant.plotPattern("vehicle", 30);

%% Footprint on the road
centre = ant.boresightGroundPoint(mountHeight);
x = linspace(0, 2 * centre(1), 301);
y = linspace(-centre(1), centre(1), 301).';
figure;
imagesc(x, y, 10 * log10(ant.groundGain(mountHeight, x, y)), max(gainDb(:)) + [-30 0]);
set(gca, "YDir", "normal");
axis equal tight;
colorbar;
xlabel("Along track x [m]");
ylabel("Cross track y [m]");
title("One-way gain on the road [dBi]");
%% Optional: table from a full-wave simulated patch antenna
% Antenna Toolbox solves the patch with the method of moments (slow). Its
% patch radiates along +z, so tilt it to point along +x (the RadarAntenna
% boresight) before sampling the pattern.
useAntennaToolbox = false;
if useAntennaToolbox
    patch = design(patchMicrostrip, fc);
    patch.Tilt = 90;
    patch.TiltAxis = [0 1 0];
    tableAz = -180:2:180;
    tableEl = -90:2:90;
    patchGainDb = pattern(patch, fc, tableAz, tableEl);
    patchAntenna = RadarAntenna(tableAz, tableEl, patchGainDb, zeros(size(patchGainDb)), ...
                                fc, depression, yaw, 0);
    fprintf("Patch boresight gain %.1f dBi\n", patchAntenna.gainDb(0, 0));
end
