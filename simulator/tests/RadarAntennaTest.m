classdef RadarAntennaTest < matlab.unittest.TestCase
    %RADARANTENNATEST Unit tests for the RadarAntenna component.
    %   Run with run_tests("RadarAntennaTest") from simulator/. Tables are
    %   generated from analytic patterns, so expected values are textbook
    %   results rather than stored snapshots.

    properties (TestParameter)
        cosineBeamwidth = {10, 30, 60, 120}
        frequency = {10.525e9, 24.125e9, 77e9}
    end

    methods (TestClassSetup)
        function addComponentsToPath(testCase)
            components = fullfile(fileparts(mfilename("fullpath")), "..", "components");
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(components));
        end
    end

    methods (Static)
        function ant = fromTable(az, el, gainDb, mount)
            % 24.125 GHz antenna with zero phase; mount = [depression yaw roll] [deg].
            ant = RadarAntenna(az, el, gainDb, zeros(size(gainDb)), 24.125e9, mount(1), mount(2), mount(3));
        end

        function ant = fromPattern(power, peakGainDb, mount)
            % Antenna from a normalised power pattern power(az, el) on a 0.25 deg grid.
            az = -180:0.25:180;
            el = (-90:0.25:90).';
            ant = RadarAntennaTest.fromTable(az, el, peakGainDb + 10 * log10(power(az, el)), mount);
        end

        function ant = gaussian(beamwidth, mount)
            % Elliptical Gaussian beam with half-power widths [az el] [deg].
            power = @(az, el) exp(-4 * log(2) * ((az / beamwidth(1)).^2 + (el / beamwidth(2)).^2));
            ant = RadarAntennaTest.fromPattern(power, 20, mount);
        end
    end

    %% Table interpolation
    methods (Test)
        function reproducesTableAtGridNodes(testCase)
            rng(1);
            az = -180:5:175;
            el = -90:5:90;
            gainDb = 20 * rand(numel(el), numel(az)) - 10;
            ant = RadarAntennaTest.fromTable(az, el, gainDb, [45 0 0]);
            testCase.verifyEqual(ant.gainDb(az, el.'), gainDb, AbsTol=1e-9);
        end

        function interpolatesBetweenNodes(testCase)
            ant = RadarAntennaTest.gaussian([8 14], [45 0 0]);
            az = -20:0.37:20;
            el = (-20:0.41:20).';
            expected = 100 * exp(-4 * log(2) * ((az / 8).^2 + (el / 14).^2));
            testCase.verifyEqual(ant.gain(az, el), expected, AbsTol=0.2);
            testCase.verifyEqual(ant.gain([-4 4 0 0], [0 0 -7 7]), 50 * ones(1, 4), RelTol=1e-3);
        end

        function carriesPhase(testCase)
            az = -180:5:180;
            el = -90:5:90;
            ant = RadarAntenna(az, el, 6 + zeros(numel(el), numel(az)), repmat(az, numel(el), 1), ...
                24.125e9, 45, 0, 0);
            testCase.verifyEqual(rad2deg(angle(ant.field(30, -10))), 30, AbsTol=1e-9);
            testCase.verifyEqual(abs(ant.field(35, 10)), 10^(6 / 20), RelTol=1e-12);
        end

        function broadcastsInputs(testCase)
            ant = RadarAntennaTest.gaussian([10 10], [45 0 0]);
            az = -40:10:40;
            el = (-20:5:20).';
            G = ant.gain(az, el);
            testCase.verifySize(G, [numel(el) numel(az)]);
            testCase.verifyEqual(G(3, 4), ant.gain(az(4), el(3)));
            testCase.verifySize(ant.gain(zeros(2, 3, 4), 0), [2 3 4]);
        end

        function noRadiationOutsideTable(testCase)
            az = -90:90;
            el = -45:45;
            ant = RadarAntennaTest.fromTable(az, el, zeros(numel(el), numel(az)), [45 0 0]);
            testCase.verifyEqual(ant.gain([0 135 -170 0], [0 0 0 60]), [1 0 0 0]);
        end

        function fullCircleTableWrapsAtSeam(testCase)
            az = -180:10:170;
            el = -90:10:90;
            ant = RadarAntennaTest.fromTable(az, el, repmat(-abs(az) / 10, numel(el), 1), [45 0 0]);
            testCase.verifyEqual(ant.gain(175, 0), ant.gain(-175, 0), RelTol=1e-12);
            testCase.verifyEqual(ant.gainDb(180, 0), -18, AbsTol=1e-12);
            testCase.verifyGreaterThan(ant.gain(175, 0), 0);
        end

        function invalidTablesThrow(testCase)
            az = -10:10;
            el = -5:5;
            gainDb = zeros(numel(el), numel(az));
            phaseDeg = zeros(size(gainDb));
            make = @(az, el, gainDb, phaseDeg, frequency, depression) ...
                RadarAntenna(az, el, gainDb, phaseDeg, frequency, depression, 0, 0);
            testCase.verifyError(@() make(az, el, gainDb.', phaseDeg, 24e9, 45), "RadarAntenna:invalidTable");
            testCase.verifyError(@() make(fliplr(az), el, gainDb, phaseDeg, 24e9, 45), "RadarAntenna:invalidTable");
            testCase.verifyError(@() make(az, el, gainDb, 0, 24e9, 45), "RadarAntenna:invalidTable");
            testCase.verifyError(@() make(az * 20, el, gainDb, phaseDeg, 24e9, 45), "MATLAB:validators:mustBeInRange");
            testCase.verifyError(@() make(az, el, NaN(size(gainDb)), phaseDeg, 24e9, 45), "MATLAB:validators:mustBeNonNan");
            testCase.verifyError(@() make(az, el, gainDb, phaseDeg, 0, 45), "MATLAB:validators:mustBePositive");
            testCase.verifyError(@() make(az, el, gainDb, phaseDeg, 24e9, 95), "MATLAB:validators:mustBeInRange");
        end

        function allArgumentsAreRequired(testCase)
            az = -10:10;
            el = -5:5;
            gainDb = zeros(numel(el), numel(az));
            args = {az, el, gainDb, zeros(size(gainDb)), 24e9, 45, 0, 0};
            for n = 0:numel(args) - 1
                testCase.verifyError(@() RadarAntenna(args{1:n}), "MATLAB:minrhs");
            end
            ant = RadarAntenna(args{:});
            testCase.verifyError(@() ant.plotPattern("vehicle"), "MATLAB:minrhs");
            testCase.verifyError(@() ant.plotPattern(), "MATLAB:minrhs");
        end

        function wavelengthMatchesFrequency(testCase, frequency)
            ant = RadarAntenna([-1 1], [-1 1], zeros(2), zeros(2), frequency, 45, 0, 0);
            testCase.verifyEqual(ant.Wavelength, 299792458 / frequency, RelTol=1e-15);
        end
    end

    %% Directivity
    methods (Test)
        function isotropicTableHasUnitDirectivity(testCase)
            ant = RadarAntennaTest.fromTable(-180:180, -90:90, zeros(181, 361), [45 0 0]);
            testCase.verifyEqual(ant.directivity(), 1, RelTol=1e-4);
        end

        function cosineDirectivityMatchesTheory(testCase, cosineBeamwidth)
            % D = 2(n+1) for a cos^n(theta) power pattern.
            n = log(0.5) / log(cosd(cosineBeamwidth / 2));
            power = @(az, el) max(cosd(az) .* cosd(el), 0).^n;
            ant = RadarAntennaTest.fromPattern(power, 0, [45 0 0]);
            testCase.verifyEqual(ant.directivity(), 2 * (n + 1), RelTol=2e-3);
        end

        function directivityIgnoresAbsoluteLevel(testCase)
            ant = RadarAntennaTest.gaussian([10 20], [45 0 0]);
            louder = RadarAntennaTest.fromTable(ant.Azimuth, ant.Elevation, ant.GainDb + 7, [45 0 0]);
            testCase.verifyEqual(louder.directivity(), ant.directivity(), RelTol=1e-12);
        end
    end

    %% Mounting on the vehicle
    methods (Test)
        function boresightFollowsMountingAngles(testCase)
            ant = RadarAntennaTest.gaussian([10 10], [30 0 0]);
            testCase.verifyEqual(ant.boresight(), [cosd(30); 0; -sind(30)], AbsTol=1e-15);
            yawed = RadarAntennaTest.gaussian([10 10], [30 90 40]);
            testCase.verifyEqual(yawed.boresight(), [0; cosd(30); -sind(30)], AbsTol=1e-15);
            R = yawed.mountingRotation();
            testCase.verifyEqual(R.' * R, eye(3), AbsTol=1e-14);
            testCase.verifyEqual(det(R), 1, AbsTol=1e-14);
        end

        function vehicleDirectionsMapToAntennaAngles(testCase)
            ant = RadarAntennaTest.gaussian([10 10], [37 -12 25]);
            [az, el] = ant.vehicleToAntenna(5 * ant.boresight());
            testCase.verifyEqual([az el], [0 0], AbsTol=1e-12);
            R = ant.mountingRotation();
            local = [cosd(10) * cosd(20); cosd(10) * sind(20); sind(10)];
            [az, el] = ant.vehicleToAntenna(R * local);
            testCase.verifyEqual([az el], [20 10], AbsTol=1e-12);
        end

        function depressionTiltsTheVerticalPlane(testCase)
            % Without roll the antenna elevation plane is the vertical plane
            % along the direction of travel.
            ant = RadarAntennaTest.gaussian([30 10], [5 0 0]);
            [az, el] = ant.vehicleToAntenna([1; 0; 0]);
            testCase.verifyEqual([az el], [0 5], AbsTol=1e-12);
            testCase.verifyEqual(ant.gain(az, el), 50, RelTol=1e-3);
        end

        function rollSwapsBeamPlanes(testCase)
            ant = RadarAntennaTest.gaussian([5 20], [0 0 90]);
            [az, el] = ant.vehicleToAntenna([cosd(10); sind(10); 0]);
            testCase.verifyEqual([az el], [0 -10], AbsTol=1e-12);
            testCase.verifyEqual(ant.gain(az, el), 50, RelTol=1e-3);
        end

        function groundGeometry(testCase)
            height = 0.6;
            ant = RadarAntennaTest.gaussian([12 12], [35 0 0]);
            p = ant.boresightGroundPoint(height);
            testCase.verifyEqual(p, [height / tand(35); 0], AbsTol=1e-12);
            [G, slant] = ant.groundGain(height, p(1), p(2));
            testCase.verifyEqual(G, 100, RelTol=1e-9);
            testCase.verifyEqual(slant, height / sind(35), RelTol=1e-12);

            G = ant.groundGain(height, linspace(0, 3, 31), linspace(-1, 1, 21).');
            testCase.verifySize(G, [21 31]);
            testCase.verifyEqual(G, flipud(G), "RelTol", 1e-9, "Footprint is symmetric across track.");
            upward = RadarAntennaTest.gaussian([12 12], [-5 0 0]);
            testCase.verifyEqual(upward.boresightGroundPoint(height), [NaN; NaN]);
        end
    end

    %% Plots
    methods (Test)
        function cutsShowPrincipalPlanes(testCase)
            original = get(groot, "DefaultFigureVisible");
            set(groot, "DefaultFigureVisible", "off");
            testCase.addTeardown(@set, groot, "DefaultFigureVisible", original);

            ant = RadarAntennaTest.gaussian([12 20], [45 0 0]);
            cutAngles = -90:0.5:90;
            fig = ant.plotCuts(cutAngles, 60);
            testCase.addTeardown(@close, fig);
            ax = findobj(fig, "Type", "axes");
            lines = findobj(ax, "Type", "line");
            azCut = findobj(lines, "DisplayName", "Azimuth plane");
            elCut = findobj(lines, "DisplayName", "Elevation plane");
            testCase.verifyEqual(azCut.XData, cutAngles);
            testCase.verifyEqual(azCut.YData, ant.gainDb(cutAngles, 0));
            testCase.verifyEqual(elCut.YData, ant.gainDb(0, cutAngles));
            testCase.verifyEqual(ax.YLim, [-40 25], AbsTol=1e-9);
            testCase.verifyError(@() ant.plotCuts(cutAngles), "MATLAB:minrhs");
        end

        function patternPeaksAlongBoresight(testCase)
            original = get(groot, "DefaultFigureVisible");
            set(groot, "DefaultFigureVisible", "off");
            testCase.addTeardown(@set, groot, "DefaultFigureVisible", original);

            ant = RadarAntennaTest.gaussian([10 20], [35 20 0]);
            frames = ["antenna", "vehicle"];
            expected = [[1; 0; 0], ant.boresight()];
            for k = 1:2
                fig = ant.plotPattern(frames(k), 30);
                testCase.addTeardown(@close, fig);
                s = findobj(fig, "Type", "surface");
                xyz = [s.XData(:).'; s.YData(:).'; s.ZData(:).'];
                [radius, i] = max(vecnorm(xyz));
                testCase.verifyEqual(radius, 30, AbsTol=1e-9);
                testCase.verifyEqual(xyz(:, i) / radius, expected(:, k), AbsTol=1e-9);
                testCase.verifyEqual(s.CData(i), 20, AbsTol=1e-9);
            end
        end

        function patternShowsVehicleAndAntennaAxes(testCase)
            original = get(groot, "DefaultFigureVisible");
            set(groot, "DefaultFigureVisible", "off");
            testCase.addTeardown(@set, groot, "DefaultFigureVisible", original);

            ant = RadarAntennaTest.gaussian([10 20], [35 20 10]);
            R = ant.mountingRotation();
            axisTip = @(fig, tag) [findobj(fig, "Tag", tag).XData(2); ...
                findobj(fig, "Tag", tag).YData(2); findobj(fig, "Tag", tag).ZData(2)] / 52;

            fig = ant.plotPattern("vehicle", 40);
            testCase.addTeardown(@close, fig);
            testCase.verifyEqual(axisTip(fig, "vehicle-x"), [1; 0; 0], AbsTol=1e-12);
            testCase.verifyEqual(axisTip(fig, "antenna-x"), ant.boresight(), AbsTol=1e-12);
            testCase.verifyEqual(axisTip(fig, "antenna-z"), R(:, 3), AbsTol=1e-12);
            testCase.verifyNumElements(findobj(fig, "Tag", "origin"), 1);

            fig = ant.plotPattern("antenna", 40);
            testCase.addTeardown(@close, fig);
            testCase.verifyEqual(axisTip(fig, "antenna-x"), [1; 0; 0], AbsTol=1e-12);
            testCase.verifyEqual(axisTip(fig, "vehicle-z"), R.' * [0; 0; 1], AbsTol=1e-12);
        end
    end
end
