close all

%%% GSS ground echo waveform model

%% Radar parameters
fc = 300e6;                       % Carrier frequency [Hz]
fs = 20*fc;                     % Sample frequency [Hz]
c = physconst('LightSpeed');    % Speed of light [m/s]
lambda = c/fc;                  % Carrier wavelength [m]
prf = 1000;                     % PRF dictates how often we update the platform/geometry phase
pri = 1/prf;
samplesperblock = round(fs*pri);

%% Generate carrier waveform and spectrum
t = (0:samplesperblock-1)' / fs; % Continuous-wave signal over one PRI
x = sin(2*pi*fc*t);              % Transmit waveform [vector]
sampleCount = numel(x);         % Number of simulated samples
carrierfft = fft(x);            
carrierspectrum = abs(fftshift(carrierfft)) / sampleCount;
frequencyAxis = (-floor(sampleCount/2):ceil(sampleCount/2)-1) * (fs / sampleCount);


%% Antenna radiation pattern setup
azang = -180:180;
elang = -90:90;

% Narrow beam along the x-direction (azimuth = 0 degrees)
n = 90; % Higher values produce a narrower beam
azpattern = sind(azang).^n;
azpattern(azpattern < 1e-10) = 1e-10;
elpattern = cosd(elang).^n;
elpattern(elpattern < 1e-10) = 1e-10;
magpattern = mag2db(elpattern(:) * azpattern);
phasepattern = zeros(size(magpattern));

antenna = phased.CustomAntennaElement('AzimuthAngles', azang, ...
                                      'ElevationAngles', elang, ...
                                      'MagnitudePattern', magpattern, ...
                                      'PhasePattern', phasepattern);

% 3D representation of the antenna radiation pattern
pattern(antenna,fc,-180:180,-90:90,'CoordinateSystem','polar','Type','powerdb')

sensor = phased.ULA('NumElements', 1,'Element', antenna);
%sensor = phased.URA('NumElements', 1,'Element', antenna);

%% Environment setup
v = 1e8; % Speed parallel to ground [m/s]
h = 1; % Antenna height above ground [m]
theta = 45; % Incidence angle from horizontal [degrees]
gammaDb = 0; % Ground gamma value at specified frequency [dB]

% ConstantGammaModel: sigma = 10^(gamma/10)*sin(theta) [theta = grazing
% angle]
gclutter = constantGammaClutter('Sensor', sensor, ...
                                'PropagationSpeed', c, ...
                                'PRF', prf, ...
                                'Gamma', gammaDb, ...
                                'EarthModel', 'Flat', ...
                                'SampleRate', fs, ...
                                'OperatingFrequency', fc, ...
                                'OutputFormat', 'Samples', ...
                                'TransmitSignalInputPort', true, ...
                                'PlatformHeight', h, ...
                                'PlatformSpeed', v, ...
                                'PlatformDirection', [90;0], ...
                                'MountingAngles', [45, 0, 0], ...
                                'NumSamples', samplesperblock);
%% Simulate the echo
yGround = gclutter(x);

maxDopplerHz = 2*v/lambda; % Maximum monostatic Doppler for radial motion

%% Plot results
% Time axis for the clutter return
Tground = (0:size(yGround,1)-1)' / fs;
echoCount = numel(yGround);
echoFFT = fft(yGround);
echoSpectrum = abs(fftshift(echoFFT)) / echoCount;
echoFrequencyAxis = (-floor(echoCount/2):ceil(echoCount/2)-1) * (fs / echoCount);

figure;
subplot(2,2,1);
plot(t, x);
grid on;
xlabel('Time [s]');
ylabel('Amplitude');
title('Transmitted waveform');
xlim([0 20/fc]);

subplot(2,2,2);
plot(frequencyAxis/1e9, carrierspectrum);
grid on;
xlabel('Frequency [GHz]');
ylabel('Magnitude');
title('Transmitted spectrum');

subplot(2,2,3);
plot(Tground, yGround);
grid on;
xlabel('Time [s]');
ylabel('Received signal amplitude');
title('Echo waveform');
xlim([0 20/fc]);

subplot(2,2,4);
plot(echoFrequencyAxis/1e6, echoSpectrum);
grid on;
xlabel('Frequency [MHz]');
ylabel('Magnitude');
title('Echo spectrum');
xlim(([fc - 1.25*maxDopplerHz, fc + 1.25*maxDopplerHz]) / 1e6);
