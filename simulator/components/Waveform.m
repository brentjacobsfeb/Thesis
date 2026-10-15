classdef Waveform
    %WAVEFORM Signal passed along the measurement chain.
    %
    %   wf = Waveform(samples, sampleRate, carrierFrequency) holds the
    %   complex baseband samples [sqrt(W)] of a signal sampled at sampleRate
    %   [Hz] around carrierFrequency [Hz]; |samples|^2 is the instantaneous
    %   power [W]. A carrierFrequency of 0 holds a signal at its own
    %   frequency, such as the output of a mixer.
    %
    %   wf = Waveform(..., Label=label) names the signal; plots use the
    %   label as their title (default "Waveform").
    %
    %   Every component takes a waveform and returns a modified copy, so
    %   each stage of the chain can be kept and inspected:
    %
    %       radiated = ant.transmit(tx);
    %       echo     = ground.scatter(radiated, height);
    %       rx       = ant.receive(echo);
    %       plot([tx, echo, rx])
    %
    %   State
    %     A waveform is "conducted" on a wire or antenna port, "radiated" by
    %     the antenna in Radiator, or "incident" on the radar after
    %     scattering. The samples of a radiated waveform are the signal at
    %     the port of the radiating antenna; those of an incident waveform
    %     are the signal an isotropic antenna would receive, and
    %     IncidentField gives the signal of any other antenna.

    properties
        Samples (:,1) double {mustBeNonempty, mustBeFinite} = 0 % Complex baseband [sqrt(W)]
        SampleRate (1,1) double {mustBePositive, mustBeFinite} = 1  % [Hz]
        CarrierFrequency (1,1) double {mustBeNonnegative, mustBeFinite}  % [Hz]
        Label (1,1) string = "Waveform"     % Name of the stage, used in plots
        Radiator = []       % Antenna radiating the waveform (empty unless radiated)
        IncidentField = []  % samples = IncidentField(antenna): port signal of a
                            % receiving antenna at the radar (empty unless incident)
    end

    properties (Dependent)
        State           % "conducted", "radiated" or "incident"
        Time            % Sample instants [s]
        Power           % Mean power [W]
        Wavelength      % Carrier wavelength [m]
    end

    methods
        function obj = Waveform(samples, sampleRate, carrierFrequency, options)
            arguments
                samples
                sampleRate
                carrierFrequency
                options.Label (1,1) string = "Waveform"
            end
            obj.Samples = samples;
            obj.SampleRate = sampleRate;
            obj.CarrierFrequency = carrierFrequency;
            obj.Label = options.Label;
        end

        function s = get.State(obj)
            if ~isempty(obj.IncidentField)
                s = "incident";
            elseif ~isempty(obj.Radiator)
                s = "radiated";
            else
                s = "conducted";
            end
        end

        function t = get.Time(obj)
            t = (0:numel(obj.Samples) - 1).' / obj.SampleRate;
        end

        function p = get.Power(obj)
            p = mean(abs(obj.Samples).^2);
        end

        function lambda = get.Wavelength(obj)
            lambda = 299792458 / obj.CarrierFrequency;
        end

        function [P, f] = spectrum(obj)
            %SPECTRUM Power per frequency bin.
            %   [P, f] = spectrum(wf) returns the power P [W] in each bin of
            %   the discrete Fourier transform and the bin frequencies f [Hz]
            %   relative to the carrier, from -SampleRate/2 upwards. The bins
            %   are SampleRate/numel(Samples) wide and sum(P) equals Power.
            n = numel(obj.Samples);
            P = fftshift(abs(fft(obj.Samples)).^2) / n^2;
            f = (-floor(n / 2):ceil(n / 2) - 1).' * obj.SampleRate / n;
        end

        function fig = plot(obj, options)
            %PLOT Time signal and spectrum of one or more waveforms.
            %   plot(wf) draws the time signal (left) and the spectrum
            %   (right) of each waveform in the array wf, one row per
            %   waveform, so plot([tx, echo, rx]) compares the stages of the
            %   chain. The options are those of plotTime and plotSpectrum.
            arguments
                obj
                options.TimeLimits (1,2) double = [0 Inf]
                options.DynamicRange (1,1) double {mustBePositive} = 80
            end
            fig = figure("Name", strjoin([obj.Label], ", "));
            layout = tiledlayout(fig, numel(obj), 2, "TileSpacing", "compact");
            for k = 1:numel(obj)
                obj(k).plotTime(Parent=nexttile(layout), TimeLimits=options.TimeLimits);
                obj(k).plotSpectrum(Parent=nexttile(layout), DynamicRange=options.DynamicRange);
            end
        end

        function [ax, unit] = plotTime(obj, options)
            %PLOTTIME In-phase and quadrature samples against time.
            %   [ax, unit] = plotTime(wf) plots in new axes and returns them
            %   with the time unit [s] of the x axis, so that xline(ax, t /
            %   unit) marks the instant t [s]. A real signal gets one line.
            %
            %   Options (name-value)
            %     Parent     - axes to plot in (default: a new figure).
            %     TimeLimits - shown interval [s] (default [0 Inf], all).
            arguments
                obj (1,1)
                options.Parent = []
                options.TimeLimits (1,2) double = [0 Inf]
            end
            ax = options.Parent;
            if isempty(ax)
                ax = axes(figure("Name", obj.Label));
            end
            span = min(options.TimeLimits(2), numel(obj.Samples) / obj.SampleRate) - options.TimeLimits(1);
            [unit, prefix] = siScale(span);
            t = obj.Time / unit;
            if any(imag(obj.Samples))
                plot(ax, t, real(obj.Samples), "DisplayName", "I");
                hold(ax, "on");
                plot(ax, t, imag(obj.Samples), "DisplayName", "Q");
                hold(ax, "off");
                legend(ax);
            else
                plot(ax, t, real(obj.Samples), "DisplayName", "Signal");
            end
            grid(ax, "on");
            xlim(ax, options.TimeLimits / unit);
            xlabel(ax, "Time [" + prefix + "s]");
            ylabel(ax, "Amplitude [\surdW]");
            title(ax, obj.Label, "Time signal");
        end

        function [ax, unit] = plotSpectrum(obj, options)
            %PLOTSPECTRUM Power per frequency bin [dBm] around the carrier.
            %   [ax, unit] = plotSpectrum(wf) plots the spectrum (see
            %   spectrum) in new axes and returns them with the frequency
            %   unit [Hz] of the x axis, so that xline(ax, f / unit) marks
            %   the frequency f [Hz] relative to the carrier.
            %
            %   Options (name-value)
            %     Parent       - axes to plot in (default: a new figure).
            %     DynamicRange - [dB] shown below the peak (default 80);
            %                    weaker bins are drawn at that floor.
            arguments
                obj (1,1)
                options.Parent = []
                options.DynamicRange (1,1) double {mustBePositive} = 80
            end
            ax = options.Parent;
            if isempty(ax)
                ax = axes(figure("Name", obj.Label));
            end
            [P, f] = obj.spectrum();
            [unit, prefix] = siScale(obj.SampleRate / 2);
            powerDbm = 10 * log10(P / 1e-3);
            peakDbm = max(powerDbm);
            % Bins below the shown range sit on its floor, so that a tone
            % between empty (-Inf dBm) bins still draws a line.
            plot(ax, f / unit, max(powerDbm, peakDbm - options.DynamicRange), "DisplayName", "Spectrum");
            grid(ax, "on");
            xlim(ax, [-1 1] * obj.SampleRate / 2 / unit);
            if isfinite(peakDbm)
                ylim(ax, peakDbm + [-options.DynamicRange 5]);
            end
            if obj.CarrierFrequency > 0
                [carrierUnit, carrierPrefix] = siScale(obj.CarrierFrequency);
                xlabel(ax, sprintf("Frequency - %g %sHz [%sHz]", ...
                    obj.CarrierFrequency / carrierUnit, carrierPrefix, prefix));
            else
                xlabel(ax, "Frequency [" + prefix + "Hz]");
            end
            ylabel(ax, "Power [dBm]");
            title(ax, obj.Label, sprintf("Spectrum: %.1f dBm in total, %.3g Hz bins", ...
                10 * log10(obj.Power / 1e-3), obj.SampleRate / numel(obj.Samples)));
        end
    end
end

function [unit, prefix] = siScale(value)
% Power of 1000 between 1e-9 and 1e9, and its SI prefix, that shows value
% with one to three digits before the decimal point.
exponent = 3 * floor(log10(abs(value)) / 3 + 1e-9);
if ~isfinite(exponent)
    exponent = 0;
end
exponent = min(max(exponent, -9), 9);
prefixes = ["n", "\mu", "m", "", "k", "M", "G"];
prefix = prefixes(exponent / 3 + 4);
unit = 10^exponent;
end
