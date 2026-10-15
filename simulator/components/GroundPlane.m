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

        function fD = dopplerShift(obj, antenna, height, x, y)
            %DOPPLERSHIFT Doppler shift [Hz] of the echo from road points (x, y) [m].
            %   The antenna sits at height [m] above (0, 0); x and y can be
            %   arrays of compatible size, as in RadarAntenna.groundGain.
            x = x + 0 * y;
            y = y + 0 * x;
            slantRange = sqrt(x.^2 + y.^2 + height^2);
            fD = -2 * obj.Velocity * x ./ (slantRange * antenna.Wavelength);
        end

        function rx = backscatter(obj, antenna, height, tx, fs, options)
            %BACKSCATTER Echo of the road received by the antenna.
            %   rx = backscatter(ground, antenna, height, tx, fs) returns the
            %   complex baseband signal [sqrt(W)] at the antenna port when the
            %   antenna, mounted at height [m] above the road, transmits the
            %   complex baseband signal tx [sqrt(W)] sampled at fs [Hz]. rx has
            %   the size of tx.
            %
            %   The road is a set of point scatterers with random positions and
            %   circular complex Gaussian amplitudes that move with the ground.
            %   The antenna transmits and receives (monostatic), so scatterer k
            %   at slant range R_k(t) contributes, by the radar equation,
            %       tx(t) * lambda * F_k^2 * sqrt(sigma_k) / ((4*pi)^1.5 * R_k^2) * exp(-j*4*pi*R_k/lambda)
            %   with F_k the complex antenna field towards it and sigma_k =
            %   sigma0 * dA its radar cross section. The changing phase
            %   4*pi*R_k/lambda produces the Doppler shift. The delay 2R/c (a
            %   few ns) is neglected in the envelope of tx, which holds for CW
            %   and slowly modulated signals.
            %
            %   Every call draws a new road from the global random stream; call
            %   rng(seed) first for a repeatable echo.
            %
            %   Options (name-value)
            %     NumScatterers - scatterers on the illuminated road at any
            %                     instant (default 2000).
            %     DynamicRange  - [dB] the road is populated where the echo
            %                     power per unit area is within DynamicRange
            %                     of its peak (default 40).
            arguments
                obj
                antenna (1,1) RadarAntenna
                height (1,1) double {mustBePositive, mustBeFinite}
                tx {mustBeVector, mustBeNumeric, mustBeFinite}
                fs (1,1) double {mustBePositive, mustBeFinite}
                options.NumScatterers (1,1) double {mustBeInteger, mustBePositive} = 2000
                options.DynamicRange (1,1) double {mustBePositive, mustBeFinite} = 40
            end
            lambda = antenna.Wavelength;
            [xLit, yLit] = GroundPlane.illuminatedArea(antenna, height, options.DynamicRange);

            % Fastest echo: the far end of the illuminated area, nearest the x axis.
            xFar = max(abs(xLit));
            yNear = max(yLit(1), 0) + max(-yLit(2), 0);
            maxDoppler = abs(obj.dopplerShift(antenna, height, xFar, yNear));
            if maxDoppler >= fs / 2
                warning("GroundPlane:aliasing", ...
                    "Doppler shifts up to %.0f Hz alias at fs = %.0f Hz; raise fs above %.0f Hz.", ...
                    maxDoppler, fs, 2 * maxDoppler);
            end

            % Scatter the road over the strip that passes through the
            % illuminated area during the signal.
            t = (0:numel(tx) - 1).' / fs;
            travel = obj.Velocity * t(end);
            x0Limits = xLit - [max(travel, 0), min(travel, 0)];
            density = options.NumScatterers / (diff(xLit) * diff(yLit));
            count = round(density * diff(x0Limits) * diff(yLit));
            x0 = x0Limits(1) + diff(x0Limits) * rand(count, 1);
            y0 = yLit(1) + diff(yLit) * rand(count, 1);
            % sigma_k = gamma * sin(grazing) / density with sin(grazing) = height / R_k.
            amplitude = sqrt(obj.Gamma * height / density) * complex(randn(count, 1), randn(count, 1)) / sqrt(2);

            % Work in time chunks over the scatterers inside the illuminated
            % area, moving at most a tenth of its length per chunk.
            chunkLength = max(1, floor(min(2^20 / options.NumScatterers, ...
                0.1 * diff(xLit) * fs / abs(obj.Velocity))));
            channel = zeros(numel(tx), 1);
            for first = 1:chunkLength:numel(tx)
                idx = first:min(first + chunkLength - 1, numel(tx));
                shift = obj.Velocity * t(idx([1 end]));
                active = x0 + max(shift) >= xLit(1) & x0 + min(shift) <= xLit(2);
                x = x0(active).' + obj.Velocity * t(idx);
                y = y0(active).' + 0 * x;
                slantRange = sqrt(x.^2 + y.^2 + height^2);
                [az, el] = antenna.vehicleToAntenna([x(:).'; y(:).'; -height * ones(1, numel(x))]);
                F = reshape(antenna.field(az, el), size(x));
                echoes = F.^2 ./ slantRange.^2.5 .* exp(-4j * pi / lambda * slantRange);
                channel(idx) = echoes * amplitude(active);
            end
            rx = lambda / (4 * pi)^1.5 * reshape(channel, size(tx)) .* tx;
        end
    end

    methods (Static, Access = private)
        function [xLimits, yLimits] = illuminatedArea(antenna, height, dynamicRange)
            % Bounding box [m] of the road where the echo power per unit area,
            % proportional to G^2 * sigma0 / R^4 ~ G^2 * sin(grazing)^5, is
            % within dynamicRange [dB] of its peak. The road is sampled every
            % 0.25 deg in vehicle azimuth and grazing angle.
            az = -180:0.25:180;
            grazing = (0.25:0.25:90).';
            groundRange = height ./ tand(grazing);
            x = groundRange .* cosd(az);
            y = groundRange .* sind(az);
            power = antenna.groundGain(height, x, y).^2 .* sind(grazing).^5;
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
