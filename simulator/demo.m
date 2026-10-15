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
patternRange = 30;      % Dynamic range of the 3D pattern [dB]
patternFig = ant.plotPattern("vehicle", patternRange);

%% Footprint on the road
centre = ant.boresightGroundPoint(mountHeight);
x = linspace(0, 2 * centre(1), 301);
y = linspace(-centre(1), centre(1), 301).';
roadGainDb = 10 * log10(ant.groundGain(mountHeight, x, y));
figure;
imagesc(x, y, roadGainDb, max(gainDb(:)) + [-30 0]);
set(gca, "YDir", "normal");
axis equal tight;
colorbar;
xlabel("Along track x [m]");
ylabel("Cross track y [m]");
title("One-way gain on the road [dBi]");

% Illuminated road in the 3D pattern. The pattern radius is in dB, so the
% road is scaled to put the boresight ground point at the tip of the main
% lobe (radius patternRange); it shares the colour scale of the pattern.
ax = findobj(patternFig, "Type", "axes");
scale = patternRange / norm([centre; mountHeight]);
hold(ax, "on");
surf(ax, scale * x, scale * y, -scale * mountHeight + 0 * roadGainDb, max(roadGainDb, ax.CLim(1)), ...
    "EdgeColor", "none", "FaceAlpha", 0.85, "Tag", "road");
text(ax, scale * centre(1), scale * centre(2), -scale * mountHeight, "  illuminated road", ...
    "VerticalAlignment", "top");
hold(ax, "off");

%% Ground echo and Doppler shift
vehicleSpeed = 20;      % Speed over the road [m/s]
gammaDb = -10;          % Ground reflectivity [dB]
txPower = 10e-3;        % Transmit power [W]
fs = 20e3;              % Baseband sample rate [Hz]
duration = 0.2;         % Observation time [s]

ground = GroundPlane(gammaDb, -vehicleSpeed);           % The road moves backwards past the vehicle
tx = sqrt(txPower) * ones(round(duration * fs), 1);     % CW carrier at complex baseband
rx = ground.backscatter(ant, mountHeight, tx, fs);

% Speed from the spectral centroid. The nearer, steeper part of the beam
% returns more power, which pulls the centroid below the boresight Doppler
% shift, so the estimate reads slightly low.
n = numel(rx);
f = (-floor(n / 2):ceil(n / 2) - 1).' * fs / n;
spectrum = fftshift(abs(fft(rx)).^2) / n^2;             % Power per bin [W]
centroid = sum(f .* spectrum) / sum(spectrum);
boresightDoppler = ground.dopplerShift(ant, mountHeight, centre(1), centre(2));
b = ant.boresight();
estimatedSpeed = centroid * ant.Wavelength / (2 * b(1));
fprintf("Echo power %.1f dBm, Doppler centroid %.0f Hz (boresight %.0f Hz), estimated speed %.2f m/s\n", ...
    10 * log10(mean(abs(rx).^2) / 1e-3), centroid, boresightDoppler, estimatedSpeed);

figure;
subplot(2, 1, 1);
t = (0:n - 1) / fs;
plot(t * 1e3, real(rx), t * 1e3, imag(rx));
grid on;
xlim([0 10]);
xlabel("Time [ms]");
ylabel("Echo [\surdW]");
legend("I", "Q");
title("Received echo");
subplot(2, 1, 2);
plot(f / 1e3, 10 * log10(spectrum / 1e-3));
xline(boresightDoppler / 1e3, "--", "Boresight Doppler");
grid on;
xlabel("Frequency [kHz]");
ylabel("Power [dBm]");
title("Echo spectrum");

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
