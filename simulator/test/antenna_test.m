close all;
clear;

%% Radar and target parameters
c = physconst('LightSpeed');
fc = 10e9;                         % Continuous-wave carrier frequency [Hz]
sampleRate = 100e3;                % Complex-baseband sample rate [Hz]
observationTime = 0.05;            % Observation duration [s]
initialRange = 1000;               % Plane range at transmit time zero [m]
planeSpeed = 30;                   % Positive speed means moving away [m/s]
reflectionCoefficient = 1;

%% Custom antenna power pattern
% Replace these angle grids and the pattern matrix with the desired pattern.
azimuthAngles = -180:180;
elevationAngles = -90:90;
azimuthBeamwidth = 30;             % Full width at half maximum [degrees]
elevationBeamwidth = 30;           % Full width at half maximum [degrees]

azimuthPower = exp(-4 * log(2) * (azimuthAngles / azimuthBeamwidth).^2);
elevationPower = exp(-4 * log(2) * (elevationAngles / elevationBeamwidth).^2);
patternPower = elevationPower(:) * azimuthPower;
patternPower = patternPower / max(patternPower(:));

% Direction from the antenna boresight to the plane.
targetAzimuth = 10;                % [degrees]
targetElevation = 0;               % [degrees]
antennaPowerGain = interp2(azimuthAngles, elevationAngles, patternPower, ...
                           targetAzimuth, targetElevation, 'linear', 0);

%% Continuous-wave monostatic echo
% Calculate each received sample using the transmit time that produced it.
t = (0:round(observationTime * sampleRate) - 1)' / sampleRate;
transmitEnvelope = ones(size(t));  % CW envelope; RF carrier is exp(j*2*pi*fc*t)

% A signal transmitted at tTx reaches the moving plane and returns at tRx.
% This gives the exact two-way Doppler factor (c-v)/(c+v), rather than the
% small-v approximation 1-2v/c.
receiveTime = t;
transmitTime = ((c - planeSpeed) .* receiveTime - 2 * initialRange) ...
               ./ (c + planeSpeed);
hasEcho = transmitTime >= 0;
scatterTime = (c .* receiveTime - initialRange) ./ (c + planeSpeed);
scatterRange = initialRange + planeSpeed .* scatterTime;

receivedEnvelope = complex(zeros(size(t)));
receivedEnvelope(hasEcho) = reflectionCoefficient * antennaPowerGain ...
    .* (initialRange ./ scatterRange(hasEcho)).^2 ...
    .* exp(1j * 2 * pi * fc ...
          .* (transmitTime(hasEcho) - receiveTime(hasEcho)));

theoreticalDopplerHz = fc * ((c - planeSpeed) / (c + planeSpeed) - 1);

%% Time- and frequency-domain plots
% The complex envelopes are coherently downconverted to baseband, so the
% received carrier offset appears directly at the Doppler frequency.
transmitSpectrum = abs(fftshift(fft(transmitEnvelope)));
receiveSpectrum = abs(fftshift(fft(receivedEnvelope)));
frequencyAxis = (-floor(numel(t) / 2):ceil(numel(t) / 2) - 1)' ...
                 * (sampleRate / numel(t));
transmitSpectrumDb = 20 * log10(transmitSpectrum / max(transmitSpectrum) + eps);
receiveSpectrumDb = 20 * log10(receiveSpectrum / max(receiveSpectrum) + eps);

figure;
subplot(2, 2, 1);
plot(t, real(transmitEnvelope), 'DisplayName', 'Transmit');
hold on;
plot(t, real(receivedEnvelope), 'DisplayName', 'Received (real)');
grid on;
xlabel('Time [s]');
ylabel('Baseband amplitude');
title('Transmit and received CW envelopes');
legend('Location', 'best');

subplot(2, 2, 2);
plot(t, real(receivedEnvelope), 'DisplayName', 'In-phase');
hold on;
plot(t, imag(receivedEnvelope), 'DisplayName', 'Quadrature');
grid on;
xlabel('Time [s]');
ylabel('Received amplitude');
title('Received complex-baseband waveform');
legend('Location', 'best');

subplot(2, 2, 3);
plot(frequencyAxis, transmitSpectrumDb, 'DisplayName', 'Transmit');
hold on;
plot(frequencyAxis, receiveSpectrumDb, 'DisplayName', 'Received');
grid on;
xlabel('Frequency relative to carrier [Hz]');
ylabel('Normalized magnitude [dB]');
title('Baseband spectrum');
xlim(5 * abs(theoreticalDopplerHz) * [-1 1]);
legend('Location', 'best');

subplot(2, 2, 4);
imagesc(azimuthAngles, elevationAngles, 10 * log10(patternPower + eps));
axis xy;
colorbar;
xlabel('Azimuth [degrees]');
ylabel('Elevation [degrees]');
title(sprintf('Custom antenna power pattern (target at %.1f deg az, %.1f deg el)', ...
              targetAzimuth, targetElevation));
caxis([-40 0]);

%% 3D radiation pattern
[azimuthGrid, elevationGrid] = meshgrid(deg2rad(azimuthAngles), ...
                                        deg2rad(elevationAngles));
[patternX, patternY, patternZ] = sph2cart(azimuthGrid, elevationGrid, patternPower);
patternPowerDb = 10 * log10(patternPower + eps);

figure;
surf(patternX, patternY, patternZ, patternPowerDb, 'EdgeColor', 'none');
hold on;
targetRadius = interp2(azimuthAngles, elevationAngles, patternPower, ...
                       targetAzimuth, targetElevation, 'linear', 0);
[targetX, targetY, targetZ] = sph2cart(deg2rad(targetAzimuth), ...
                                       deg2rad(targetElevation), targetRadius);
plot3(targetX, targetY, targetZ, 'r.', 'MarkerSize', 24);
axis equal;
grid on;
xlabel('X (boresight) [normalized]');
ylabel('Y [normalized]');
zlabel('Z [normalized]');
title('3D normalized antenna power pattern');
colormap turbo;
colorbar;
caxis([-40 0]);
view(35, 25);

fprintf('Theoretical monostatic Doppler shift: %.3f Hz\n', theoreticalDopplerHz);
fprintf('Antenna power gain toward plane: %.3f (%.2f dB)\n', ...
        antennaPowerGain, 10 * log10(antennaPowerGain + eps));