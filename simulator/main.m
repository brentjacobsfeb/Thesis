close all

%%% GSS ground echo waveform model

%% Radar parameters
fc = 24e9;                      % Carrier frequency [Hz]
fs = 10*fc;                    % Sample frequency [Hz]
c = physconst('LightSpeed');    % Speed of light [m/s]
lambda = c/fc;                  % Wavelength [m]
h = 1;                          % Antenna height above ground [m]
theta = pi/4;                   % Incidence angle from horizontal [rad]
R = h/sin(theta);               % Direct range from ground [m]
tau = 2*R/c;                    % Round-trip delay [s]

%% Generate carrier waveform
t = (0:1/fs:10/fc);             % Time points
x = sin(2*pi*fc*t).';                    % transmit waveform used by the clutter model
sampleCount = numel(x);
frequencyFFT = fft(x);
frequencySpectrum = abs(fftshift(frequencyFFT)) / sampleCount;
frequencyAxis = (-floor(sampleCount/2):ceil(sampleCount/2)-1) * (fs / sampleCount);


%% Antenna setup
% Theoretical antenna beam: broad in azimuth, narrow in elevation,
% with peak gain along the boresight direction. This is enough to
% populate a valid CustomAntennaElement for the clutter simulation.
fc = 1e9;
azang = -180:180;
elang = -90:90;
magpattern = mag2db(repmat(cosd(elang)',1,numel(azang)));
phasepattern = zeros(size(magpattern));

antenna = phased.CustomAntennaElement('AzimuthAngles',azang, ...
    'ElevationAngles',elang,'MagnitudePattern',magpattern, ...
    'PhasePattern',phasepattern);
resp = antenna(fc,[20;30]);
pattern(antenna,fc,20,-90:90,'CoordinateSystem','polar','Type','powerdb')


sensor = phased.ULA('NumElements', 1,'Element', antenna);

%% Environment setup
Rmin = 50;                 % Minimum clutter range [m]
Rmax = 100;                % Maximum clutter range [m]
gammaDb = -20;             % Ground reflectivity gamma value [dB]
PtERP = 1;                 % ERP [W]

%% Configure radar clutter model
% Capture the pulse plus the expected round-trip delay to Rmax.
Nrecord = numel(x) + ceil(2*Rmax/c*fs);

gc = constantGammaClutter( ...
    'Sensor', sensor, ...
    'OperatingFrequency', fc, ...
    'SampleRate', fs, ...
    'PlatformHeight', h, ...
    'PlatformSpeed', 0, ...
    'Gamma', gammaDb, ...
    'TransmitSignalInputPort', true, ...
    'TransmitERP', PtERP, ...
    'OutputFormat', 'Samples', ...
    'NumSamples', Nrecord, ...
    'ClutterMinRange', Rmin, ...
    'ClutterMaxRange', Rmax);

gc.MountingAngles = [0, -theta*180/pi, 0];

%% Simulate the received ground return
yGround = gc(x);

%% Plot results
% Time axis for the clutter return
Tground = (0:size(yGround,1)-1)' / fs;

figure;
subplot(2,2,1);
plot(t, x);
grid on;
xlabel('Time [s]');
ylabel('Amplitude');
title('Transmitted waveform');

subplot(2,2,2);
plot(frequencyAxis/1e9, frequencySpectrum);
grid on;
xlabel('Frequency [GHz]');
ylabel('Magnitude');
title('Transmitted spectrum');

subplot(2,2,3);
plot(Tground, yGround);
grid on;
xlabel('Time [s]');
ylabel('Received magnitude');
title('Echo waveform');


subplot(2,2,4);
echoCount = numel(yGround);
echoFFT = fft(yGround);
echoSpectrum = abs(fftshift(echoFFT)) / echoCount;
echoFrequencyAxis = (-floor(echoCount/2):ceil(echoCount/2)-1) * (fs / echoCount);

plot(echoFrequencyAxis/1e9, echoSpectrum);
grid on;
xlabel('Frequency [GHz]');
ylabel('Magnitude');
title('Echo spectrum');
