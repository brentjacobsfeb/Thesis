classdef RadarAntenna
    %RADARANTENNA Tabulated antenna pattern mounted on the vehicle.
    %
    %   ant = RadarAntenna(az, el, gainDb, phaseDeg, frequency, depression, yaw, roll)
    %   interpolates the gain table gainDb [dBi] of size numel(el) x numel(az),
    %   given on the azimuth grid az [deg] and elevation grid el [deg] of the
    %   antenna frame. Directions outside the table do not radiate.
    %
    %     phaseDeg   - phase pattern [deg], same size as gainDb.
    %     frequency  - operating frequency [Hz].
    %     depression - boresight angle below the horizon [deg].
    %     yaw        - boresight azimuth from the direction of travel [deg].
    %     roll       - rotation about the boresight [deg].
    %
    %   Coordinate frames
    %     Antenna frame: x = boresight, y = left, z = up. Azimuth is measured
    %       from x towards y, elevation from the x-y plane towards z (the
    %       Phased Array System Toolbox convention).
    %     Vehicle frame: x = direction of travel, y = left, z = up. The
    %       antenna frame follows by turning the boresight by Yaw about z,
    %       pitching it down by Depression and rolling it by Roll.

    properties (SetAccess = private)
        Azimuth         % Table azimuth grid [deg]
        Elevation       % Table elevation grid [deg]
        GainDb          % Table gain [dBi]
        PhaseDeg        % Table phase [deg]
        Frequency       % [Hz]
        Depression      % [deg]
        Yaw             % [deg]
        Roll            % [deg]
    end

    properties (Dependent)
        Wavelength      % [m]
    end

    properties (Access = private)
        Interpolant
    end

    methods
        function obj = RadarAntenna(az, el, gainDb, phaseDeg, frequency, depression, yaw, roll)
            arguments
                az {mustBeVector, mustBeInRange(az, -180, 180)}
                el {mustBeVector, mustBeInRange(el, -90, 90)}
                gainDb {mustBeNumeric, mustBeNonNan}
                phaseDeg {mustBeNumeric, mustBeFinite}
                frequency (1,1) double {mustBePositive, mustBeFinite}
                depression (1,1) double {mustBeInRange(depression, -90, 90)}
                yaw (1,1) double {mustBeFinite}
                roll (1,1) double {mustBeFinite}
            end
            az = reshape(az, 1, []);
            el = reshape(el, 1, []);
            if numel(az) < 2 || numel(el) < 2 || any(diff(az) <= 0) || any(diff(el) <= 0)
                error("RadarAntenna:invalidTable", "The angle grids need at least two strictly increasing values.");
            end
            if ~isequal(size(gainDb), [numel(el) numel(az)]) || ~isequal(size(phaseDeg), size(gainDb))
                error("RadarAntenna:invalidTable", "gainDb and phaseDeg must be numel(el) x numel(az).");
            end
            obj.Azimuth = az;
            obj.Elevation = el;
            obj.GainDb = gainDb;
            obj.PhaseDeg = phaseDeg;
            obj.Frequency = frequency;
            obj.Depression = depression;
            obj.Yaw = yaw;
            obj.Roll = roll;

            F = 10.^(gainDb / 20) .* exp(1j * deg2rad(phaseDeg));
            % Wrap the azimuth seam when the table covers the full circle.
            seamGap = 360 - (az(end) - az(1));
            if seamGap > 0 && seamGap <= max(diff(az))
                az = [az(end) - 360, az, az(1) + 360];
                F = [F(:, end), F, F(:, 1)];
            end
            obj.Interpolant = griddedInterpolant({el, az}, F, "linear", "none");
        end

        function lambda = get.Wavelength(obj)
            lambda = 299792458 / obj.Frequency;
        end

        function F = field(obj, az, el)
            %FIELD Complex far field sqrt(gain)*exp(j*phase) in the antenna frame.
            %   az, el [deg] can be any arrays of compatible size.
            el = el + 0 * az;
            az = az + 0 * el;
            F = reshape(obj.Interpolant(el(:), mod(az(:) + 180, 360) - 180), size(az));
            F(isnan(F)) = 0;
        end

        function G = gain(obj, az, el)
            %GAIN Linear power gain in the antenna frame (az, el in degrees).
            G = abs(obj.field(az, el)).^2;
        end

        function G = gainDb(obj, az, el)
            %GAINDB Power gain [dBi] in the antenna frame (az, el in degrees).
            G = 10 * log10(obj.gain(az, el));
        end

        function D = directivity(obj)
            %DIRECTIVITY Peak directivity (linear) implied by the shape of the table.
            %   Compare 10*log10(D) with the peak of GainDb to check that a
            %   table is physically consistent; the difference is the loss.
            az = -180:0.25:180;
            el = (-90:0.25:90).';
            P = obj.gain(az, el);
            radiated = trapz(deg2rad(az), trapz(deg2rad(el), P .* cosd(el), 1));
            D = 4 * pi * max(P(:)) / radiated;
        end

        function fig = plotCuts(obj, cutAngles, dynamicRange)
            %PLOTCUTS Gain along the principal planes through boresight.
            %   plotCuts(ant, cutAngles, dynamicRange) plots the gain [dBi] in
            %   the azimuth (el = 0) and elevation (az = 0) planes at
            %   cutAngles [deg], showing dynamicRange [dB] below the peak.
            arguments
                obj
                cutAngles {mustBeVector, mustBeFinite}
                dynamicRange (1,1) double {mustBePositive}
            end
            peakDb = max(obj.GainDb(:));
            fig = figure("Name", "Antenna pattern cuts");
            ax = axes(fig);
            plot(ax, cutAngles, obj.gainDb(cutAngles, 0), cutAngles, obj.gainDb(0, cutAngles));
            grid(ax, "on");
            ylim(ax, peakDb + [-dynamicRange 5]);
            xlabel(ax, "Angle from boresight [deg]");
            ylabel(ax, "Gain [dBi]");
            legend(ax, "Azimuth plane", "Elevation plane");
        end

        function fig = plotPattern(obj, frame, dynamicRange)
            %PLOTPATTERN 3D view of the radiation pattern.
            %   plotPattern(ant, frame, dynamicRange) draws the gain as a
            %   surface whose radius is the gain [dB] above the peak minus
            %   dynamicRange [dB]. frame is "antenna" or "vehicle"; the
            %   vehicle frame shows the pattern as mounted (x = direction of
            %   travel, z = up).
            arguments
                obj
                frame (1,1) string {mustBeMember(frame, ["antenna", "vehicle"])}
                dynamicRange (1,1) double {mustBePositive}
            end
            az = -180:180;
            el = (-90:90).';
            G = obj.gainDb(az, el);
            floorDb = max(G(:)) - dynamicRange;
            r = max(G - floorDb, 0);
            X = r .* cosd(el) .* cosd(az);
            Y = r .* cosd(el) .* sind(az);
            Z = r .* sind(el) + 0 * az;
            if frame == "vehicle"
                xyz = obj.mountingRotation() * [X(:).'; Y(:).'; Z(:).'];
                X = reshape(xyz(1, :), size(X));
                Y = reshape(xyz(2, :), size(Y));
                Z = reshape(xyz(3, :), size(Z));
            end

            fig = figure("Name", "Antenna radiation pattern");
            ax = axes(fig);
            surf(ax, X, Y, Z, max(G, floorDb), "EdgeColor", "none");
            axis(ax, "equal");
            clim(ax, [floorDb, max(G(:))]);
            cb = colorbar(ax);
            cb.Label.String = "Gain [dBi]";
            xlabel(ax, "x");
            ylabel(ax, "y");
            zlabel(ax, "z");
            title(ax, sprintf("Radiation pattern (%s frame), radius = gain above %.1f dBi", ...
                frame, floorDb));

            % Vehicle (solid) and antenna (dashed) axes; both frames share the origin.
            R = obj.mountingRotation();
            if frame == "vehicle"
                frameAxes = {eye(3), R};
            else
                frameAxes = {R.', eye(3)};
            end
            frameNames = ["vehicle", "antenna"];
            lineStyles = ["-", "--"];
            axisNames = ["x", "y", "z"];
            axisColours = [0.85 0.15 0.15; 0.15 0.65 0.15; 0.15 0.3 0.9];
            axisLength = 1.3 * dynamicRange;
            hold(ax, "on");
            for k = 1:3
                labels = axisNames(k) + "_{" + frameNames + "}";
                if norm(frameAxes{1}(:, k) - frameAxes{2}(:, k)) < 1e-9
                    labels = [strjoin(labels, " = "), ""];    % coinciding axes share one label
                end
                for f = 1:2
                    tip = axisLength * frameAxes{f}(:, k);
                    plot3(ax, [0 tip(1)], [0 tip(2)], [0 tip(3)], lineStyles(f), ...
                        "Color", axisColours(k, :), "LineWidth", 2, "Tag", frameNames(f) + "-" + axisNames(k));
                    text(ax, tip(1), tip(2), tip(3), labels(f), "Color", axisColours(k, :), "FontWeight", "bold");
                end
            end
            plot3(ax, 0, 0, 0, "ko", "MarkerFaceColor", "k", "Tag", "origin");
            text(ax, 0, 0, 0, "  origin", "VerticalAlignment", "top");
            hold(ax, "off");
            view(ax, 3);
            rotate3d(fig, "on");
        end

        function R = mountingRotation(obj)
            %MOUNTINGROTATION Columns are the antenna x, y, z axes in the vehicle frame.
            cy = cosd(obj.Yaw); sy = sind(obj.Yaw);
            cp = cosd(obj.Depression); sp = sind(obj.Depression);
            cr = cosd(obj.Roll); sr = sind(obj.Roll);
            R = [cy -sy 0; sy cy 0; 0 0 1] * [cp 0 sp; 0 1 0; -sp 0 cp] * [1 0 0; 0 cr -sr; 0 sr cr];
        end

        function b = boresight(obj)
            %BORESIGHT Unit boresight vector in the vehicle frame.
            R = obj.mountingRotation();
            b = R(:, 1);
        end

        function [az, el] = vehicleToAntenna(obj, dirs)
            %VEHICLETOANTENNA Antenna-frame angles [deg] of vehicle-frame directions (3xN).
            local = obj.mountingRotation().' * dirs;
            az = atan2d(local(2, :), local(1, :));
            el = atan2d(local(3, :), vecnorm(local(1:2, :)));
        end

        function [G, slantRange] = groundGain(obj, height, x, y)
            %GROUNDGAIN Linear gain towards points (x, y) [m] on a flat road.
            %   The antenna sits at height [m] above (0, 0); x is along the
            %   direction of travel, y to the left. x and y can be arrays of
            %   compatible size.
            x = x + 0 * y;
            y = y + 0 * x;
            [az, el] = obj.vehicleToAntenna([x(:).'; y(:).'; -height * ones(1, numel(x))]);
            G = reshape(obj.gain(az, el), size(x));
            slantRange = sqrt(x.^2 + y.^2 + height^2);
        end

        function p = boresightGroundPoint(obj, height)
            %BORESIGHTGROUNDPOINT Road point [x; y] hit by the boresight (NaN if none).
            b = obj.boresight();
            if b(3) >= 0
                p = [NaN; NaN];
            else
                p = height / -b(3) * b(1:2);
            end
        end
    end
end
