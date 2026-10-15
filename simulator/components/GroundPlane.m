classdef GroundPlane
    %GROUNDPLANE Rough dielectric road that scatters the radar signal back to the antenna.
    %
    %   ground = GroundPlane(permittivity, rmsHeight, correlationLength, velocity)
    %   models the road as a randomly rough surface of a dielectric, moving
    %   past the antenna.
    %
    %     permittivity      - relative complex permittivity eps' - j*eps''
    %                         of the road at the carrier frequency.
    %     rmsHeight         - standard deviation s of the surface height [m].
    %     correlationLength - correlation length l of the surface height [m].
    %     velocity          - ground velocity along the vehicle x axis [m/s].
    %                         A vehicle driving forward at speed v sees
    %                         velocity = -v.
    %
    %   ground = GroundPlane(..., Name=Value) also sets
    %     Polarization        - "VV" (default) or "HH": the antenna
    %                           polarization, vertical or horizontal to the
    %                           plane of incidence.
    %     CorrelationFunction - "exponential" (default) or "gaussian": the
    %                           shape exp(-r/l) or exp(-r^2/l^2) of the
    %                           surface height correlation.
    %
    %   In the measurement chain, scatter turns a Waveform radiated by the
    %   antenna into the echo incident back on it.
    %
    %   Scattering
    %     The surface height is a Gaussian random process. The road returns
    %     the wave in two ways:
    %       Specular - the mean (coherent) field reflects as from a smooth
    %                  road, with the Fresnel reflection coefficient (see
    %                  fresnel) weakened by the roughness (see
    %                  specularReflection). Only the road straight below
    %                  the antenna reflects it back to the radar.
    %       Diffuse  - the roughness scatters the rest in all directions,
    %                  with the backscatter coefficient sigma0 of the
    %                  integral equation model (see sigma0). This is the
    %                  Doppler-shifted echo of the road under the beam.
    %     The polarization is taken relative to the plane of incidence of
    %     every road point, which holds for a narrow beam; the antenna model
    %     itself carries no polarization.
    %
    %   Geometry
    %     The road is the plane z = -height of the vehicle frame (x = direction
    %     of travel, y = left, z = up, antenna at the origin; see
    %     RadarAntenna). A road point at slant range R returns its echo with
    %     the Doppler shift fD = -2 * (dR/dt) / wavelength, so the road ahead
    %     of a vehicle driving forward gives a positive shift.

    properties (SetAccess = private)
        Permittivity        % Relative complex permittivity eps' - j*eps''
        RmsHeight           % Surface height standard deviation [m]
        CorrelationLength   % Surface height correlation length [m]
        CorrelationFunction % "exponential" or "gaussian"
        Polarization        % "VV" or "HH"
        Velocity            % Ground velocity along vehicle x [m/s]
    end

    methods
        function obj = GroundPlane(permittivity, rmsHeight, correlationLength, velocity, options)
            arguments
                permittivity (1,1) double {mustBeFinite}
                rmsHeight (1,1) double {mustBePositive, mustBeFinite}
                correlationLength (1,1) double {mustBePositive, mustBeFinite}
                velocity (1,1) double {mustBeReal, mustBeFinite}
                options.Polarization (1,1) string {mustBeMember(options.Polarization, ["VV", "HH"])} = "VV"
                options.CorrelationFunction (1,1) string ...
                    {mustBeMember(options.CorrelationFunction, ["exponential", "gaussian"])} = "exponential"
            end
            if real(permittivity) <= 1 || imag(permittivity) > 0
                error("GroundPlane:invalidPermittivity", ...
                    "The permittivity needs a real part above 1 and a negative imaginary part (eps' - j*eps'').");
            end
            obj.Permittivity = permittivity;
            obj.RmsHeight = rmsHeight;
            obj.CorrelationLength = correlationLength;
            obj.CorrelationFunction = options.CorrelationFunction;
            obj.Polarization = options.Polarization;
            obj.Velocity = velocity;
        end

        function R = fresnel(obj, grazing)
            %FRESNEL Fresnel reflection coefficient of the smooth road at grazing angles [deg].
            %   For the polarization of the ground, with theta = 90 deg -
            %   grazing the incidence angle and eps the permittivity,
            %       VV: R = (eps*cos(theta) - q) / (eps*cos(theta) + q)
            %       HH: R = (cos(theta) - q) / (cos(theta) + q)
            %   with q = sqrt(eps - sin(theta)^2). VV vanishes at the
            %   Brewster angle; both tend to -1 at grazing incidence.
            arguments
                obj
                grazing {mustBeNumeric, mustBeReal, mustBeInRange(grazing, 0, 90)}
            end
            cosTheta = sind(grazing);
            q = sqrt(obj.Permittivity - cosd(grazing).^2);
            if obj.Polarization == "VV"
                R = (obj.Permittivity * cosTheta - q) ./ (obj.Permittivity * cosTheta + q);
            else
                R = (cosTheta - q) ./ (cosTheta + q);
            end
        end

        function r = specularReflection(obj, source)
            %SPECULARREFLECTION Coherent reflection coefficient of the rough road at normal incidence.
            %   r = specularReflection(ground, source) is the Fresnel
            %   coefficient (1 - sqrt(eps)) / (1 + sqrt(eps)) of the smooth
            %   road times the roughness loss exp(-2*(k*s)^2) of a Gaussian
            %   height distribution, with k = 2*pi/wavelength of source, the
            %   RadarAntenna or Waveform. At normal incidence both
            %   polarizations reflect alike.
            k = 2 * pi / source.Wavelength;
            sqrtEps = sqrt(obj.Permittivity);
            r = (1 - sqrtEps) / (1 + sqrtEps) * exp(-2 * (k * obj.RmsHeight)^2);
        end

        function s0 = sigma0(obj, source, grazing)
            %SIGMA0 Diffuse backscatter coefficient [m^2/m^2] at grazing angles [deg].
            %   s0 = sigma0(ground, source, grazing) uses the wavenumber
            %   k = 2*pi/wavelength of source, the RadarAntenna or Waveform;
            %   grazing can be any array.
            %
            %   Single-scattering integral equation model (IEM; Fung, Li and
            %   Chen, IEEE Trans. Geosci. Remote Sens. 30(2), 1992). With
            %   theta = 90 deg - grazing the incidence angle and
            %   kz = k*cos(theta),
            %       sigma0 = k^2/2 * exp(-2*kz^2*s^2) * sum_{n>=1} s^(2n) * |I_n|^2 * W_n(2*k*sin(theta)) / n!
            %       I_n    = (2*kz)^n * f * exp(-kz^2*s^2) + kz^n * F/2
            %   where the Kirchhoff and complementary field coefficients
            %   follow from the Fresnel coefficient R (see fresnel),
            %       VV: f =  2*R/cos(theta)
            %           F =  2*sin(theta)^2 * (1+R)^2 / cos(theta) * (1 - 1/eps + (eps - sin(theta)^2 - eps*cos(theta)^2) / (eps^2*cos(theta)^2))
            %       HH: f = -2*R/cos(theta)
            %           F = -2*sin(theta)^2 * (1+R)^2 / cos(theta) * (eps - 1) / cos(theta)^2
            %   and W_n is the Fourier transform, divided by 2*pi, of the
            %   n-th power of the correlation function:
            %       exponential: W_n(K) = (l/n)^2 * (1 + (K*l/n)^2)^(-3/2)
            %       Gaussian:    W_n(K) = l^2/(2*n) * exp(-(K*l)^2/(4*n))
            %   For k*s << 1 the first term is the small perturbation model;
            %   for large k*s and k*l with the Gaussian correlation the sum
            %   tends to geometric optics. The model is accurate up to about
            %   k*s = 3.
            arguments
                obj
                source
                grazing {mustBeNumeric, mustBeReal, mustBeInRange(grazing, 0, 90)}
            end
            k = 2 * pi / source.Wavelength;
            s = obj.RmsHeight;
            l = obj.CorrelationLength;
            er = obj.Permittivity;
            s0 = zeros(size(grazing));
            lit = grazing > 0;      % Nothing scatters back at grazing incidence
            cosTheta = sind(grazing(lit));
            sinTheta2 = cosd(grazing(lit)).^2;
            R = obj.fresnel(grazing(lit));
            if obj.Polarization == "VV"
                f = 2 * R ./ cosTheta;
                F = 2 * sinTheta2 .* (1 + R).^2 ./ cosTheta ...
                    .* (1 - 1 / er + (er - sinTheta2 - er * cosTheta.^2) ./ (er^2 * cosTheta.^2));
            else
                f = -2 * R ./ cosTheta;
                F = -2 * sinTheta2 .* (1 + R).^2 * (er - 1) ./ cosTheta.^3;
            end
            kzs = k * s * cosTheta;
            K = 2 * k * sqrt(sinTheta2);

            % The weights (2*kz*s)^(2n) * exp(-4*kz^2*s^2) / n! of the
            % Kirchhoff term are Poisson with mean 4*kz^2*s^2, so the terms
            % stop contributing past n = (2*k*s)^2. Each term is computed
            % in log scale, where its factors cannot overflow.
            total = zeros(size(cosTheta));
            for n = 1:1000
                logFactorial = gammaln(n + 1) / 2;
                a = exp(n * log(2 * kzs) - 2 * kzs.^2 - logFactorial);
                b = exp(n * log(kzs) - kzs.^2 - logFactorial);
                if obj.CorrelationFunction == "exponential"
                    W = (l / n)^2 * (1 + (K * l / n).^2).^-1.5;
                else
                    W = l^2 / (2 * n) * exp(-(K * l).^2 / (4 * n));
                end
                term = abs(a .* f + b .* F / 2).^2 .* W;
                total = total + term;
                if n > (2 * k * s)^2 && all(term <= 1e-12 * total)
                    break
                end
            end
            s0(lit) = k^2 / 2 * total;
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
            %   With tx the radiated samples, lambda the carrier wavelength of
            %   wf, and Ft and Fr the complex fields of the transmitting and
            %   receiving antennas towards a direction (Fr = 1 for the
            %   isotropic antenna), the port signal of a receiving antenna at
            %   the radar is the sum of
            %     Specular - the mirror image of the antenna in the road,
            %                2*height below it, which returns
            %                    tx(t) * lambda * Ft * Fr * r / (4*pi * 2*height) * exp(-j*4*pi*height/lambda)
            %                with Ft, Fr straight down and r the coherent
            %                reflection coefficient (see specularReflection).
            %                Its path does not change as the road moves:
            %                it has no Doppler shift.
            %     Diffuse  - a set of point scatterers with random positions
            %                and circular complex Gaussian amplitudes that
            %                move with the ground. Scatterer k at slant
            %                range R_k(t) returns, by the radar equation,
            %                    tx(t) * lambda * Ft_k * Fr_k * sqrt(sigma_k) / ((4*pi)^1.5 * R_k^2) * exp(-j*4*pi*R_k/lambda)
            %                with sigma_k = sigma0(grazing_k) * dA its radar
            %                cross section (see sigma0). The changing phase
            %                4*pi*R_k/lambda produces the Doppler shift.
            %   The delay 2R/c (a few ns) is neglected in the envelope of tx,
            %   which holds for CW and slowly modulated signals.
            %
            %   Every call draws a new road from the global random stream; call
            %   rng(seed) first for a repeatable echo. Every antenna that
            %   receives the returned waveform sees the same road.
            %
            %   Options (name-value)
            %     NumScatterers - scatterers on the illuminated road at any
            %                     instant (default 2000).
            %     DynamicRange  - [dB] the road is populated where the
            %                     diffuse echo power per unit area at the
            %                     isotropic antenna is within DynamicRange
            %                     of its peak (default 40).
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
            lambda = wf.Wavelength;
            ks = 2 * pi / lambda * obj.RmsHeight;
            if ks > 3
                warning("GroundPlane:outsideValidity", ...
                    "k*s = %.1f: the IEM backscatter is accurate up to about k*s = 3.", ks);
            end
            fs = wf.SampleRate;
            [xLit, yLit] = obj.illuminatedArea(wf, height, options.DynamicRange);

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
            % sigma_k = sigma0(grazing_k) / density; sqrt(sigma0) follows the
            % scatterer from a table in grazing angle.
            road.Amplitude = sqrt(1 / density) * complex(randn(count, 1), randn(count, 1)) / sqrt(2);
            grazing = (0:0.01:90).';
            road.SqrtSigma0 = griddedInterpolant(grazing, sqrt(obj.sigma0(wf, grazing)), "linear", "nearest");
            road.Specular = obj.specularReflection(wf) * lambda / (8 * pi * height) * exp(-4j * pi * height / lambda);
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

    methods (Access = private)
        function [xLimits, yLimits] = illuminatedArea(obj, radiated, height, dynamicRange)
            % Bounding box [m] of the road where the diffuse echo power per
            % unit area at an isotropic receiving antenna, proportional to
            % G * sigma0 / R^4 ~ G * sigma0 * sin(grazing)^4, is within
            % dynamicRange [dB] of its peak. The road is sampled every
            % 0.25 deg in vehicle azimuth and grazing angle.
            az = -180:0.25:180;
            grazing = (0.25:0.25:90).';
            groundRange = height ./ tand(grazing);
            x = groundRange .* cosd(az);
            y = groundRange .* sind(az);
            power = radiated.Radiator.groundGain(height, x, y) .* obj.sigma0(radiated, grazing) .* sind(grazing).^4;
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
% when empty), for the echo of the radiated Waveform from road (see
% GroundPlane.scatter): the specular reflection straight below the antenna
% plus the scatterers, in time chunks over those inside the illuminated
% area.
lambda = radiated.Wavelength;
t = radiated.Time;
down = [0; 0; -1];
specular = road.Specular * fieldTowards(radiated.Radiator, down);
if ~isempty(receiver)
    specular = specular * fieldTowards(receiver, down);
end
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
    sqrtSigma0 = road.SqrtSigma0(asind(road.Height ./ slantRange(:)));
    echoes = reshape(F(:) .* sqrtSigma0, size(x)) ./ slantRange.^2 .* exp(-4j * pi / lambda * slantRange);
    channel(idx) = echoes * road.Amplitude(active);
end
samples = (lambda / (4 * pi)^1.5 * channel + specular) .* radiated.Samples;
end

function F = fieldTowards(antenna, directions)
% Complex field of antenna towards vehicle-frame directions (3xN).
[az, el] = antenna.vehicleToAntenna(directions);
F = antenna.field(az, el);
end
