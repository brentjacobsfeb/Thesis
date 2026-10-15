classdef GroundPlane
    %GROUNDPLANE Flat road that scatters the radar signal back to the antenna.
    %
    %   ground = GroundPlane(gammaDb, velocity) models the road as a
    %   distributed scatterer with the constant-gamma backscatter coefficient
    %   sigma0 = gamma * sin(grazing angle), moving past the antenna.
    %
    %     gammaDb  - ground reflectivity gamma [dB].
    %     velocity - ground velocity along the vehicle x axis [m/s]. A vehicle
    %                driving forward at speed v sees velocity = -v.
    %
    %   In the measurement chain, scatter turns a Waveform radiated by the
    %   antenna into the echo incident back on it.
    %
    %   Geometry
    %     The road is the plane z = -height of the vehicle frame (x = direction
    %     of travel, y = left, z = up, antenna at the origin; see
    %     RadarAntenna). A road point at slant range R returns its echo with
    %     the Doppler shift fD = -2 * (dR/dt) / wavelength, so the road ahead
    %     of a vehicle driving forward gives a positive shift.

    properties (SetAccess = private)
        GammaDb         % Reflectivity [dB]
        Velocity        % Ground velocity along vehicle x [m/s]
    end

    properties (Dependent)
        Gamma           % Reflectivity (linear)
    end

    methods
        function obj = GroundPlane(gammaDb, velocity)
            arguments
                gammaDb (1,1) double {mustBeReal, mustBeFinite}
                velocity (1,1) double {mustBeReal, mustBeFinite}
            end
            obj.GammaDb = gammaDb;
            obj.Velocity = velocity;
        end

        function g = get.Gamma(obj)
            g = 10^(obj.GammaDb / 10);
        end

        function s0 = sigma0(obj, grazing)
            %SIGMA0 Backscatter coefficient [m^2/m^2] at grazing angles [deg].
            s0 = obj.Gamma * sind(grazing);
        end

        function fD = dopplerShift(obj, source, height, x, y)
            %DOPPLERSHIFT Doppler shift [Hz] of the echo from road points (x, y) [m].
            %   source is the RadarAntenna or Waveform whose wavelength sets
            %   the shift. The antenna sits at height [m] above (0, 0); x and
            %   y can be arrays of compatible size, as in
            %   RadarAntenna.groundGain.
            x = x + 0 * y;
            y = y + 0 * x;
            slantRange = sqrt(x.^2 + y.^2 + height^2);
            fD = -2 * obj.Velocity * x ./ (slantRange * source.Wavelength);
        end

        function wf = scatter(obj, wf, height, options)
            %SCATTER Echo of the road arriving back at the radar.
            %   wf = scatter(ground, wf, height) returns the incident Waveform
            %   that reaches the radar when the Waveform wf, radiated by an
            %   antenna mounted at height [m] above the road (see
            %   RadarAntenna.transmit), hits the road. Its samples are the
            %   signal [sqrt(W)] an isotropic antenna at the radar would
            %   receive; RadarAntenna.receive gives the signal of a real one.
            %
            %   The road is a set of point scatterers with random positions and
            %   circular complex Gaussian amplitudes that move with the ground.
            %   Scatterer k at slant range R_k(t) contributes, by the radar
            %   equation,
            %       tx(t) * lambda * Ft_k * Fr_k * sqrt(sigma_k) / ((4*pi)^1.5 * R_k^2) * exp(-j*4*pi*R_k/lambda)
            %   to the port signal of a receiving antenna at the radar, with
            %   tx the radiated samples, Ft_k and Fr_k the complex fields of
            %   the transmitting and receiving antennas towards the scatterer
            %   (Fr_k = 1 for the isotropic antenna) and sigma_k = sigma0 * dA
            %   its radar cross section. The changing phase 4*pi*R_k/lambda,
            %   with lambda the carrier wavelength of wf, produces the Doppler
            %   shift. The delay 2R/c (a few ns) is neglected in the envelope
            %   of tx, which holds for CW and slowly modulated signals.
            %
            %   Every call draws a new road from the global random stream; call
            %   rng(seed) first for a repeatable echo. Every antenna that
            %   receives the returned waveform sees the same road.
            %
            %   Options (name-value)
            %     NumScatterers - scatterers on the illuminated road at any
            %                     instant (default 2000).
            %     DynamicRange  - [dB] the road is populated where the echo
            %                     power per unit area at the isotropic antenna
            %                     is within DynamicRange of its peak
            %                     (default 40).
            arguments
                obj
                wf (1,1) Waveform
                height (1,1) double {mustBePositive, mustBeFinite}
                options.NumScatterers (1,1) double {mustBeInteger, mustBePositive} = 2000
                options.DynamicRange (1,1) double {mustBePositive, mustBeFinite} = 40
            end
            if wf.State ~= "radiated"
                error("GroundPlane:notRadiated", ...
                    "scatter needs a waveform radiated by an antenna, not a %s one.", wf.State);
            end
            fs = wf.SampleRate;
            [xLit, yLit] = GroundPlane.illuminatedArea(wf.Radiator, height, options.DynamicRange);

            % Fastest echo: the far end of the illuminated area, nearest the x axis.
            xFar = max(abs(xLit));
            yNear = max(yLit(1), 0) + max(-yLit(2), 0);
            maxDoppler = abs(obj.dopplerShift(wf, height, xFar, yNear));
            if maxDoppler >= fs / 2
                warning("GroundPlane:aliasing", ...
                    "Doppler shifts up to %.0f Hz alias at fs = %.0f Hz; raise fs above %.0f Hz.", ...
                    maxDoppler, fs, 2 * maxDoppler);
            end

            % Scatter the road over the strip that passes through the
            % illuminated area during the signal.
            travel = obj.Velocity * (numel(wf.Samples) - 1) / fs;
            x0Limits = xLit - [max(travel, 0), min(travel, 0)];
            density = options.NumScatterers / (diff(xLit) * diff(yLit));
            count = round(density * diff(x0Limits) * diff(yLit));
            road.X0 = x0Limits(1) + diff(x0Limits) * rand(count, 1);
            road.Y0 = yLit(1) + diff(yLit) * rand(count, 1);
            % sigma_k = gamma * sin(grazing) / density with sin(grazing) = height / R_k.
            road.Amplitude = sqrt(obj.Gamma * height / density) * complex(randn(count, 1), randn(count, 1)) / sqrt(2);
            road.XLimits = xLit;
            road.Height = height;
            road.Velocity = obj.Velocity;
            % Work in time chunks over the scatterers inside the illuminated
            % area, moving at most a tenth of its length per chunk.
            road.ChunkLength = max(1, floor(min(2^20 / options.NumScatterers, ...
                0.1 * diff(xLit) * fs / abs(obj.Velocity))));

            radiated = wf;
            wf.Samples = roadEcho(road, radiated, []);
            wf.IncidentField = @(receiver) roadEcho(road, radiated, receiver);
            wf.Radiator = [];
            wf.Label = "Road echo (isotropic antenna)";
        end
    end

    methods (Static, Access = private)
        function [xLimits, yLimits] = illuminatedArea(antenna, height, dynamicRange)
            % Bounding box [m] of the road where the echo power per unit area
            % at an isotropic receiving antenna, proportional to
            % G * sigma0 / R^4 ~ G * sin(grazing)^5, is within dynamicRange
            % [dB] of its peak. The road is sampled every 0.25 deg in vehicle
            % azimuth and grazing angle.
            az = -180:0.25:180;
            grazing = (0.25:0.25:90).';
            groundRange = height ./ tand(grazing);
            x = groundRange .* cosd(az);
            y = groundRange .* sind(az);
            power = antenna.groundGain(height, x, y) .* sind(grazing).^5;
            if ~any(power(:) > 0)
                error("GroundPlane:noIllumination", "The antenna does not illuminate the road.");
            end
            lit = power >= max(power(:)) * 10^(-dynamicRange / 10);
            % Grow by one sample so that a beam narrower than the grid keeps an area.
            lit(2:end, :) = lit(2:end, :) | lit(1:end - 1, :);
            lit(1:end - 1, :) = lit(1:end - 1, :) | lit(2:end, :);
            lit(:, 2:end) = lit(:, 2:end) | lit(:, 1:end - 1);
            lit(:, 1:end - 1) = lit(:, 1:end - 1) | lit(:, 2:end);
            xLimits = [min(x(lit)), max(x(lit))];
            yLimits = [min(y(lit)), max(y(lit))];
        end
    end
end

function samples = roadEcho(road, radiated, receiver)
% Signal at the port of receiver, a RadarAntenna at the radar (isotropic
% when empty), for the echo of the radiated Waveform from the scatterers of
% road (see GroundPlane.scatter). Works in time chunks over the scatterers
% inside the illuminated area.
lambda = radiated.Wavelength;
t = radiated.Time;
channel = zeros(numel(t), 1);
for first = 1:road.ChunkLength:numel(t)
    idx = first:min(first + road.ChunkLength - 1, numel(t));
    shift = road.Velocity * t(idx([1 end]));
    active = road.X0 + max(shift) >= road.XLimits(1) & road.X0 + min(shift) <= road.XLimits(2);
    x = road.X0(active).' + road.Velocity * t(idx);
    y = road.Y0(active).' + 0 * x;
    slantRange = sqrt(x.^2 + y.^2 + road.Height^2);
    directions = [x(:).'; y(:).'; -road.Height * ones(1, numel(x))];
    F = fieldTowards(radiated.Radiator, directions);
    if ~isempty(receiver)
        F = F .* fieldTowards(receiver, directions);
    end
    echoes = reshape(F, size(x)) ./ slantRange.^2.5 .* exp(-4j * pi / lambda * slantRange);
    channel(idx) = echoes * road.Amplitude(active);
end
samples = lambda / (4 * pi)^1.5 * channel .* radiated.Samples;
end

function F = fieldTowards(antenna, directions)
% Complex field of antenna towards vehicle-frame directions (3xN).
[az, el] = antenna.vehicleToAntenna(directions);
F = antenna.field(az, el);
end
