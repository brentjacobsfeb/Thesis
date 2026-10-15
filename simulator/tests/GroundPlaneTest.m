classdef GroundPlaneTest < matlab.unittest.TestCase
    %GROUNDPLANETEST Unit tests for the GroundPlane component.
    %   Run with run_tests("GroundPlaneTest") from simulator/. Expected
    %   values come from the radar equation and the Doppler formula
    %   integrated over the road; random tests fix the seed.

    properties (TestParameter)
        lookYaw = struct("forward", 0, "backward", 180)
        velocity = {-25, 12}
    end

    methods (TestClassSetup)
        function addComponentsToPath(testCase)
            components = fullfile(fileparts(mfilename("fullpath")), "..", "components");
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(components));
        end
    end

    methods (Static)
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

        function [power, meanDoppler] = expectedEcho(ground, ant, height, txPower)
            % Radar equation integrated over the road on a 2.5 mm grid around
            % the boresight point: mean echo power [W] and power-weighted
            % mean Doppler shift [Hz].
            p = ant.boresightGroundPoint(height);
            x = p(1) + (-1.5:2.5e-3:1.5);
            y = (-1:2.5e-3:1).';
            [G, slantRange] = ant.groundGain(height, x, y);
            w = G.^2 .* ground.sigma0(asind(height ./ slantRange)) ./ slantRange.^4;
            power = txPower * ant.Wavelength^2 / (4 * pi)^3 * sum(w, "all") * 2.5e-3^2;
            meanDoppler = sum(w .* ground.dopplerShift(ant, height, x, y), "all") / sum(w, "all");
        end
    end

    %% Reflectivity
    methods (Test)
        function constantGammaBackscatter(testCase)
            ground = GroundPlane(-10, 0);
            testCase.verifyEqual(ground.Gamma, 0.1, RelTol=1e-15);
            testCase.verifyEqual(ground.sigma0([90 30 0]), [0.1 0.05 0], AbsTol=1e-15);
        end

        function echoScalesWithGamma(testCase)
            ant = GroundPlaneTest.gaussian([12 20], [45 0 0]);
            weak = GroundPlane(-20, -10);
            strong = GroundPlane(-5, -10);
            rng(3);
            weakEcho = weak.backscatter(ant, 0.5, ones(64, 1), 20e3);
            rng(3);
            strongEcho = strong.backscatter(ant, 0.5, ones(64, 1), 20e3);
            testCase.verifyEqual(strongEcho, weakEcho * 10^(15 / 20), RelTol=1e-12);
        end

        function powerFollowsRadarEquation(testCase)
            % Mean over independent roads; 400 realisations give a 5% spread.
            ant = GroundPlaneTest.gaussian([12 20], [45 0 0]);
            ground = GroundPlane(-10, 0);
            txPower = 0.01;
            rng(4);
            rx = zeros(400, 1);
            for k = 1:numel(rx)
                rx(k) = ground.backscatter(ant, 0.5, sqrt(txPower), 1e3, NumScatterers=200);
            end
            expected = GroundPlaneTest.expectedEcho(ground, ant, 0.5, txPower);
            testCase.verifyEqual(mean(abs(rx).^2), expected, RelTol=0.15);
        end
    end

    %% Echo signal
    methods (Test)
        function staticGroundGivesConstantEcho(testCase)
            ant = GroundPlaneTest.gaussian([12 20], [45 0 0]);
            ground = GroundPlane(-10, 0);
            rx = ground.backscatter(ant, 0.5, ones(100, 1), 1e3);
            testCase.verifyNotEqual(rx(1), 0);
            testCase.verifyEqual(rx, rx(1) * ones(100, 1), RelTol=1e-12);
        end

        function echoIsLinearInTransmitSignal(testCase)
            ant = GroundPlaneTest.gaussian([12 20], [45 0 0]);
            ground = GroundPlane(-10, -10);
            n = 0:255;
            tx = (1 + n / 256) .* exp(2j * pi * 300 * n / 8e3);
            rng(6);
            cw = ground.backscatter(ant, 0.5, ones(256, 1), 8e3);
            rng(6);
            rx = ground.backscatter(ant, 0.5, 3 * tx, 8e3);
            testCase.verifySize(rx, size(tx));
            testCase.verifyEqual(rx, 3 * tx .* cw.', RelTol=1e-12);
        end

        function phasePatternActsTwoWay(testCase)
            ground = GroundPlane(-10, -10);
            rng(7);
            rx = ground.backscatter(GroundPlaneTest.gaussian([12 20], [45 0 0]), 0.5, ones(32, 1), 8e3);
            rng(7);
            shifted = ground.backscatter(GroundPlaneTest.gaussian([12 20], [45 0 0], 25), 0.5, ones(32, 1), 8e3);
            testCase.verifyEqual(shifted, rx * exp(2j * deg2rad(25)), RelTol=1e-12);
        end
    end

    %% Doppler
    methods (Test)
        function dopplerFollowsRadialVelocity(testCase, lookYaw, velocity)
            ant = GroundPlaneTest.gaussian([6 6], [40 lookYaw 0]);
            ground = GroundPlane(-10, velocity);
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
            ground = GroundPlane(-10, velocity);
            fs = 16e3;
            n = 4096;
            rng(5);
            rx = ground.backscatter(ant, height, ones(n, 1), fs);
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
            ground = GroundPlane(-10, -25);
            testCase.verifyWarning(@() ground.backscatter(ant, 0.5, ones(16, 1), 4e3), "GroundPlane:aliasing");
            testCase.verifyWarningFree(@() ground.backscatter(ant, 0.5, ones(16, 1), 20e3));
        end
    end

    %% Arguments
    methods (Test)
        function allArgumentsAreRequired(testCase)
            for n = 0:1
                args = {-10, 5};
                testCase.verifyError(@() GroundPlane(args{1:n}), "MATLAB:minrhs");
            end
            ground = GroundPlane(-10, 5);
            ant = GroundPlaneTest.gaussian([12 20], [45 0 0]);
            testCase.verifyError(@() ground.backscatter(ant, 0.5, ones(4, 1)), "MATLAB:minrhs");
        end

        function invalidArgumentsThrow(testCase)
            testCase.verifyError(@() GroundPlane(NaN, 0), "MATLAB:validators:mustBeFinite");
            testCase.verifyError(@() GroundPlane(-10, Inf), "MATLAB:validators:mustBeFinite");
            testCase.verifyError(@() GroundPlane(-10, [1 2]), "MATLAB:validation:IncompatibleSize");
            ground = GroundPlane(-10, 5);
            ant = GroundPlaneTest.gaussian([12 20], [45 0 0]);
            testCase.verifyError(@() ground.backscatter(ant, -0.5, ones(4, 1), 1e3), "MATLAB:validators:mustBePositive");
            testCase.verifyError(@() ground.backscatter(ant, 0.5, ones(4), 1e3), "MATLAB:validators:mustBeVector");
            testCase.verifyError(@() ground.backscatter("antenna", 0.5, ones(4, 1), 1e3), "MATLAB:validation:UnableToConvert");
        end

        function antennaMustSeeTheRoad(testCase)
            % A beam pointing 45 deg up with no pattern outside +-10 deg.
            ant = RadarAntenna(-10:10, -10:10, zeros(21), zeros(21), 24.125e9, -45, 0, 0);
            ground = GroundPlane(-10, 5);
            testCase.verifyError(@() ground.backscatter(ant, 0.5, ones(4, 1), 1e3), "GroundPlane:noIllumination");
        end
    end
end
