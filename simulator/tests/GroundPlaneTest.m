classdef GroundPlaneTest < matlab.unittest.TestCase
    %GROUNDPLANETEST Unit tests for the GroundPlane component.
    %   Run with run_tests("GroundPlaneTest") from simulator/. Expected
    %   values come from the Fresnel equations, the small perturbation and
    %   geometric optics limits of the scattering model, and the radar
    %   equation and Doppler formula integrated over the road; random tests
    %   fix the seed.

    properties (TestParameter)
        lookYaw = struct("forward", 0, "backward", 180)
        velocity = {-25, 12}
        polarization = struct("VV", "VV", "HH", "HH")
        correlationFunction = struct("exponential", "exponential", "gaussian", "gaussian")
    end

    methods (TestClassSetup)
        function addComponentsToPath(testCase)
            components = fullfile(fileparts(mfilename("fullpath")), "..", "components");
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(components));
        end
    end

    methods (Static)
        function ground = road(velocity, varargin)
            % Asphalt-like road: eps = 4 - 0.2j, 1 mm RMS height, 5 mm
            % correlation length; varargin are options of GroundPlane.
            ground = GroundPlane(4 - 0.2j, 1e-3, 5e-3, velocity, varargin{:});
        end

        function wf = carrier()
            % Waveform on the 24.125 GHz carrier, for its wavelength.
            wf = Waveform(0, 1, 24.125e9);
        end

        function ant = gaussian(beamwidth, mount, phaseDeg)
            % 24.125 GHz, 20 dBi Gaussian beam with half-power widths [az el]
            % [deg]; mount = [depression yaw roll] [deg].
            arguments
                beamwidth
                mount
                phaseDeg = 0
            end
            az = -180:0.25:180;
            el = (-90:0.25:90).';
            gainDb = 20 - 10 * log10(2) * 4 * ((az / beamwidth(1)).^2 + (el / beamwidth(2)).^2);
            ant = RadarAntenna(az, el, gainDb, phaseDeg + zeros(size(gainDb)), 24.125e9, mount(1), mount(2), mount(3));
        end

        function echo = radiateAndScatter(ground, ant, height, tx, fs, varargin)
            % Incident echo of tx [sqrt(W)], sampled at fs [Hz] and radiated
            % by ant at height [m]; varargin are options of scatter.
            radiated = ant.transmit(Waveform(tx, fs, ant.Frequency));
            echo = ground.scatter(radiated, height, varargin{:});
        end

        function rx = roundTrip(ground, ant, height, tx, fs, varargin)
            % Samples received by ant for the echo of tx (see radiateAndScatter).
            received = ant.receive(GroundPlaneTest.radiateAndScatter(ground, ant, height, tx, fs, varargin{:}));
            rx = received.Samples;
        end

        function [power, meanDoppler] = expectedEcho(ground, ant, height, txPower)
            % Radar equation integrated over the road on a 2.5 mm grid around
            % the boresight point: mean echo power [W] and power-weighted
            % mean Doppler shift [Hz].
            p = ant.boresightGroundPoint(height);
            x = p(1) + (-1.5:2.5e-3:1.5);
            y = (-1:2.5e-3:1).';
            [G, slantRange] = ant.groundGain(height, x, y);
            w = G.^2 .* ground.sigma0(ant, asind(height ./ slantRange)) ./ slantRange.^4;
            power = txPower * ant.Wavelength^2 / (4 * pi)^3 * sum(w, "all") * 2.5e-3^2;
            meanDoppler = sum(w .* ground.dopplerShift(ant, height, x, y), "all") / sum(w, "all");
        end
    end

    %% Surface model
    methods (Test)
        function fresnelCoefficients(testCase)
            vv = GroundPlane(4, 1e-3, 5e-3, 0, Polarization="VV");
            hh = GroundPlane(4, 1e-3, 5e-3, 0, Polarization="HH");
            testCase.verifyEqual([vv.fresnel(90) hh.fresnel(90)], [1 -1] / 3, "AbsTol", 1e-15, ...
                "Normal incidence gives (1 - sqrt(eps)) / (1 + sqrt(eps)), VV with the opposite sign convention.");
            testCase.verifyEqual([vv.fresnel(0) hh.fresnel(0)], [-1 -1], "AbsTol", 1e-15, ...
                "Grazing incidence reflects totally.");
            testCase.verifyEqual(vv.fresnel(atand(1 / 2)), 0, "AbsTol", 1e-12, ...
                "VV vanishes at the Brewster angle, atan(sqrt(eps)) from the vertical.");
            lossy = GroundPlaneTest.road(0, Polarization="HH");
            testCase.verifyLessThan(abs(lossy.fresnel(0.5:0.5:90)), 1);
        end

        function smallRoughnessMatchesSmallPerturbation(testCase, polarization, correlationFunction)
            % For k*s << 1 the backscatter is the first-order Bragg term
            % 8*k^4*s^2*cos(theta)^4*|alpha|^2*W(2*k*sin(theta)); the
            % higher-order terms add a relative O((k*s)^2).
            carrier = GroundPlaneTest.carrier();
            k = 2 * pi / carrier.Wavelength;
            er = 4 - 0.2j;
            s = 1e-6;
            l = 5e-3;
            ground = GroundPlane(er, s, l, 0, Polarization=polarization, CorrelationFunction=correlationFunction);
            theta = [0 10 30 45 60 80];
            q = sqrt(er - sind(theta).^2);
            if polarization == "VV"
                alpha = (er - 1) * (sind(theta).^2 - er * (1 + sind(theta).^2)) ./ (er * cosd(theta) + q).^2;
            else
                alpha = (cosd(theta) - q) ./ (cosd(theta) + q);
            end
            K = 2 * k * sind(theta);
            if correlationFunction == "exponential"
                W = l^2 * (1 + (K * l).^2).^-1.5;
            else
                W = l^2 / 2 * exp(-(K * l).^2 / 4);
            end
            spm = 8 * k^4 * s^2 * cosd(theta).^4 .* abs(alpha).^2 .* W;
            testCase.verifyEqual(ground.sigma0(carrier, 90 - theta), spm, RelTol=1e-3);
        end

        function largeRoughnessMatchesGeometricOptics(testCase, polarization)
            % For large k*s and k*l the Gaussian surface reflects off facets:
            % |R|^2 * exp(-tan(theta)^2 / (2*m^2)) / (2*m^2*cos(theta)^4)
            % with the RMS slope m = sqrt(2)*s/l.
            carrier = GroundPlaneTest.carrier();
            s = 10e-3;
            l = 60e-3;
            ground = GroundPlane(4 - 0.2j, s, l, 0, Polarization=polarization, CorrelationFunction="gaussian");
            theta = [0 10 20 30];
            m2 = 2 * s^2 / l^2;
            go = abs(ground.fresnel(90 - theta)).^2 .* exp(-tand(theta).^2 / (2 * m2)) ./ (2 * m2 * cosd(theta).^4);
            testCase.verifyEqual(ground.sigma0(carrier, 90 - theta), go, RelTol=0.02);
        end

        function polarizationsMeetAtNormalIncidence(testCase)
            carrier = GroundPlaneTest.carrier();
            vv = GroundPlaneTest.road(0, Polarization="VV");
            hh = GroundPlaneTest.road(0, Polarization="HH");
            testCase.verifyEqual(vv.sigma0(carrier, 90), hh.sigma0(carrier, 90), RelTol=1e-12);
            testCase.verifyEqual([vv.sigma0(carrier, 0) hh.sigma0(carrier, 0)], [0 0], ...
                "Nothing scatters back at grazing incidence.");
            testCase.verifyGreaterThan(vv.sigma0(carrier, 30:15:75), hh.sigma0(carrier, 30:15:75), ...
                "Bragg scattering from a dielectric is stronger in VV.");
        end

        function roughnessTradesSpecularForDiffuse(testCase)
            carrier = GroundPlaneTest.carrier();
            k = 2 * pi / carrier.Wavelength;
            smooth = GroundPlane(4, 0.1e-3, 5e-3, 0);
            rough = GroundPlane(4, 1e-3, 5e-3, 0);
            testCase.verifyEqual(rough.specularReflection(carrier), -1 / 3 * exp(-2 * (k * 1e-3)^2), RelTol=1e-12);
            testCase.verifyLessThan(abs(rough.specularReflection(carrier)), abs(smooth.specularReflection(carrier)));
            testCase.verifyGreaterThan(rough.sigma0(carrier, 15:15:90), smooth.sigma0(carrier, 15:15:90));
        end

        function smoothRoadReflectsLikeAMirror(testCase)
            % Looking straight down at a nearly smooth road, the echo is the
            % image of the antenna 2*height below it, with no Doppler shift
            % although the road moves.
            ant = GroundPlaneTest.gaussian([12 20], [90 0 0]);
            ground = GroundPlane(4 - 0.2j, 1e-6, 5e-3, -10);
            height = 0.5;
            lambda = ant.Wavelength;
            rng(9);
            rx = GroundPlaneTest.roundTrip(ground, ant, height, ones(64, 1), 20e3);
            mirror = 10^(20 / 10) * ground.specularReflection(ant) * lambda / (4 * pi * 2 * height) ...
                * exp(-2j * pi * 2 * height / lambda);
            testCase.verifyEqual(rx, mirror * ones(64, 1), RelTol=1e-3);
        end

        function powerFollowsRadarEquation(testCase)
            % Mean over independent roads; 400 realisations give a 5% spread.
            ant = GroundPlaneTest.gaussian([12 20], [45 0 0]);
            ground = GroundPlaneTest.road(0);
            txPower = 0.01;
            rng(4);
            rx = zeros(400, 1);
            for k = 1:numel(rx)
                rx(k) = GroundPlaneTest.roundTrip(ground, ant, 0.5, sqrt(txPower), 1e3, NumScatterers=200);
            end
            expected = GroundPlaneTest.expectedEcho(ground, ant, 0.5, txPower);
            testCase.verifyEqual(mean(abs(rx).^2), expected, RelTol=0.15);
        end
    end

    %% Echo signal
    methods (Test)
        function staticGroundGivesConstantEcho(testCase)
            ant = GroundPlaneTest.gaussian([12 20], [45 0 0]);
            ground = GroundPlaneTest.road(0);
            rx = GroundPlaneTest.roundTrip(ground, ant, 0.5, ones(100, 1), 1e3);
            testCase.verifyNotEqual(rx(1), 0);
            testCase.verifyEqual(rx, rx(1) * ones(100, 1), RelTol=1e-12);
        end

        function echoIsLinearInTransmitSignal(testCase)
            ant = GroundPlaneTest.gaussian([12 20], [45 0 0]);
            ground = GroundPlaneTest.road(-10);
            n = (0:255).';
            tx = (1 + n / 256) .* exp(2j * pi * 300 * n / 8e3);
            rng(6);
            cw = GroundPlaneTest.roundTrip(ground, ant, 0.5, ones(256, 1), 8e3);
            rng(6);
            rx = GroundPlaneTest.roundTrip(ground, ant, 0.5, 3 * tx, 8e3);
            testCase.verifySize(rx, size(tx));
            testCase.verifyEqual(rx, 3 * tx .* cw, RelTol=1e-12);
        end

        function phasePatternActsTwoWay(testCase)
            ground = GroundPlaneTest.road(-10);
            rng(7);
            rx = GroundPlaneTest.roundTrip(ground, GroundPlaneTest.gaussian([12 20], [45 0 0]), 0.5, ones(32, 1), 8e3);
            rng(7);
            shifted = GroundPlaneTest.roundTrip(ground, GroundPlaneTest.gaussian([12 20], [45 0 0], 25), 0.5, ones(32, 1), 8e3);
            testCase.verifyEqual(shifted, rx * exp(2j * deg2rad(25)), RelTol=1e-12);
        end

        function transmitAndReceivePatternsAreSeparate(testCase)
            % The incident echo carries the transmit pattern once; each
            % receiving antenna adds its own pattern.
            ground = GroundPlaneTest.road(-10);
            plain = GroundPlaneTest.gaussian([12 20], [45 0 0]);
            shifted = GroundPlaneTest.gaussian([12 20], [45 0 0], 25);
            rng(7);
            echo = GroundPlaneTest.radiateAndScatter(ground, plain, 0.5, ones(32, 1), 8e3);
            rng(7);
            shiftedEcho = GroundPlaneTest.radiateAndScatter(ground, shifted, 0.5, ones(32, 1), 8e3);
            testCase.verifyEqual(shiftedEcho.Samples, echo.Samples * exp(1j * deg2rad(25)), RelTol=1e-12);

            rx = plain.receive(echo).Samples;
            testCase.verifyEqual(shifted.receive(echo).Samples, rx * exp(1j * deg2rad(25)), RelTol=1e-12);
            testCase.verifyEqual(plain.receive(shiftedEcho).Samples, rx * exp(1j * deg2rad(25)), RelTol=1e-12);
            testCase.verifyEqual(plain.receive(echo).Samples, rx, "Receiving again sees the same road.");
        end

        function isotropicAntennaReceivesIncidentEcho(testCase)
            ground = GroundPlaneTest.road(-10);
            ant = GroundPlaneTest.gaussian([12 20], [45 0 0]);
            isotropic = RadarAntenna(-180:10:180, -90:10:90, zeros(19, 37), zeros(19, 37), 24.125e9, 45, 0, 0);
            rng(8);
            echo = GroundPlaneTest.radiateAndScatter(ground, ant, 0.5, ones(64, 1), 8e3);
            testCase.verifyEqual(isotropic.receive(echo).Samples, echo.Samples, RelTol=1e-12);
        end

        function waveformPassesThroughTheChain(testCase)
            ant = GroundPlaneTest.gaussian([12 20], [45 0 0]);
            ground = GroundPlaneTest.road(-10);
            tx = Waveform(ones(64, 1), 8e3, 24.125e9);
            radiated = ant.transmit(tx);
            testCase.verifyEqual(radiated.State, "radiated");
            testCase.verifyEqual(radiated.Samples, tx.Samples);
            echo = ground.scatter(radiated, 0.5);
            testCase.verifyEqual(echo.State, "incident");
            testCase.verifyEmpty(echo.Radiator);
            rx = ant.receive(echo);
            testCase.verifyEqual(rx.State, "conducted");
            testCase.verifyEqual([rx.SampleRate rx.CarrierFrequency], [8e3 24.125e9]);
            testCase.verifySize(rx.Samples, [64 1]);
            testCase.verifyNumElements(unique([tx.Label radiated.Label echo.Label rx.Label]), 4, ...
                "Every stage has its own label.");
        end

        function scatterNeedsRadiatedWaveform(testCase)
            ant = GroundPlaneTest.gaussian([12 20], [45 0 0]);
            ground = GroundPlaneTest.road(-10);
            tx = Waveform(ones(16, 1), 8e3, 24.125e9);
            testCase.verifyError(@() ground.scatter(tx, 0.5), "GroundPlane:notRadiated");
            echo = ground.scatter(ant.transmit(tx), 0.5);
            testCase.verifyError(@() ground.scatter(echo, 0.5), "GroundPlane:notRadiated");
        end
    end

    %% Doppler
    methods (Test)
        function dopplerFollowsRadialVelocity(testCase, lookYaw, velocity)
            ant = GroundPlaneTest.gaussian([6 6], [40 lookYaw 0]);
            ground = GroundPlaneTest.road(velocity);
            p = ant.boresightGroundPoint(0.5);
            expected = -2 * velocity * cosd(lookYaw) * cosd(40) / ant.Wavelength;
            testCase.verifyEqual(ground.dopplerShift(ant, 0.5, p(1), p(2)), expected, RelTol=1e-12);
            testCase.verifyEqual(ground.dopplerShift(ant, 0.5, 0, [-1 0 1]), [0 0 0], ...
                "Ground straight across or below the antenna has no Doppler shift.");
            fD = ground.dopplerShift(ant, 0.5, linspace(-2, 2, 5), linspace(-1, 1, 3).');
            testCase.verifySize(fD, [3 5]);
            testCase.verifyEqual(fD, flipud(fD), "Doppler shift is symmetric across track.");
        end

        function spectrumCentredOnDoppler(testCase, lookYaw, velocity)
            % The power-weighted mean Doppler shift sits slightly closer to
            % zero than the boresight value: the steeper, nearer part of the
            % beam returns more power.
            height = 0.5;
            ant = GroundPlaneTest.gaussian([6 6], [45 lookYaw 0]);
            ground = GroundPlaneTest.road(velocity);
            fs = 16e3;
            n = 4096;
            rng(5);
            rx = GroundPlaneTest.roundTrip(ground, ant, height, ones(n, 1), fs);
            spectrum = abs(fft(rx)).^2;
            f = mod((0:n - 1).' * fs / n + fs / 2, fs) - fs / 2;
            centroid = sum(f .* spectrum) / sum(spectrum);

            [~, expected] = GroundPlaneTest.expectedEcho(ground, ant, height, 1);
            testCase.verifyEqual(centroid, expected, RelTol=0.02);
            p = ant.boresightGroundPoint(height);
            boresight = ground.dopplerShift(ant, height, p(1), p(2));
            testCase.verifyEqual(sign(centroid), sign(boresight));
            testCase.verifyLessThan(abs(expected), abs(boresight));
        end

        function warnsWhenDopplerAliases(testCase)
            ant = GroundPlaneTest.gaussian([12 20], [45 0 0]);
            ground = GroundPlaneTest.road(-25);
            testCase.verifyWarning(@() GroundPlaneTest.radiateAndScatter(ground, ant, 0.5, ones(16, 1), 4e3), ...
                "GroundPlane:aliasing");
            testCase.verifyWarningFree(@() GroundPlaneTest.radiateAndScatter(ground, ant, 0.5, ones(16, 1), 20e3));
        end
    end

    %% Validity
    methods (Test)
        function warnsWhenTooRough(testCase)
            ant = GroundPlaneTest.gaussian([12 20], [45 0 0]);
            gravel = GroundPlane(4 - 0.2j, 10e-3, 20e-3, 0);
            testCase.verifyWarning(@() GroundPlaneTest.radiateAndScatter(gravel, ant, 0.5, ones(16, 1), 20e3), ...
                "GroundPlane:outsideValidity");
        end
    end

    %% Arguments
    methods (Test)
        function allArgumentsAreRequired(testCase)
            args = {4 - 0.2j, 1e-3, 5e-3, 5};
            for n = 0:3
                testCase.verifyError(@() GroundPlane(args{1:n}), "MATLAB:minrhs");
            end
            ground = GroundPlaneTest.road(5);
            ant = GroundPlaneTest.gaussian([12 20], [45 0 0]);
            radiated = ant.transmit(Waveform(ones(4, 1), 1e3, ant.Frequency));
            testCase.verifyError(@() ground.scatter(radiated), "MATLAB:minrhs");
        end

        function invalidArgumentsThrow(testCase)
            testCase.verifyError(@() GroundPlane(NaN, 1e-3, 5e-3, 0), "MATLAB:validators:mustBeFinite");
            testCase.verifyError(@() GroundPlane(0.5, 1e-3, 5e-3, 0), "GroundPlane:invalidPermittivity");
            testCase.verifyError(@() GroundPlane(4 + 0.2j, 1e-3, 5e-3, 0), "GroundPlane:invalidPermittivity");
            testCase.verifyError(@() GroundPlane(4, 0, 5e-3, 0), "MATLAB:validators:mustBePositive");
            testCase.verifyError(@() GroundPlane(4, 1e-3, -5e-3, 0), "MATLAB:validators:mustBePositive");
            testCase.verifyError(@() GroundPlane(4, 1e-3, 5e-3, Inf), "MATLAB:validators:mustBeFinite");
            testCase.verifyError(@() GroundPlane(4, 1e-3, 5e-3, [1 2]), "MATLAB:validation:IncompatibleSize");
            testCase.verifyError(@() GroundPlane(4, 1e-3, 5e-3, 0, Polarization="VH"), "MATLAB:validators:mustBeMember");
            testCase.verifyError(@() GroundPlane(4, 1e-3, 5e-3, 0, CorrelationFunction="cosine"), ...
                "MATLAB:validators:mustBeMember");
            ground = GroundPlaneTest.road(5);
            testCase.verifyError(@() ground.sigma0(GroundPlaneTest.carrier(), 95), "MATLAB:validators:mustBeInRange");
            ant = GroundPlaneTest.gaussian([12 20], [45 0 0]);
            radiated = ant.transmit(Waveform(ones(4, 1), 1e3, ant.Frequency));
            testCase.verifyError(@() ground.scatter(radiated, -0.5), "MATLAB:validators:mustBePositive");
            testCase.verifyError(@() ground.scatter(ones(4, 1), 0.5), "MATLAB:validation:UnableToConvert");
        end

        function antennaMustSeeTheRoad(testCase)
            % A beam pointing 45 deg up with no pattern outside +-10 deg.
            ant = RadarAntenna(-10:10, -10:10, zeros(21), zeros(21), 24.125e9, -45, 0, 0);
            ground = GroundPlaneTest.road(5);
            testCase.verifyError(@() GroundPlaneTest.radiateAndScatter(ground, ant, 0.5, ones(4, 1), 1e3), ...
                "GroundPlane:noIllumination");
        end
    end
end
