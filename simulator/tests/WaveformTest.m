classdef WaveformTest < matlab.unittest.TestCase
    %WAVEFORMTEST Unit tests for the Waveform passed along the chain.
    %   Run with run_tests("WaveformTest") from simulator/. Spectra are
    %   checked on tones that fall exactly on a frequency bin.

    properties (TestParameter)
        numSamples = {64, 65}
    end

    methods (TestClassSetup)
        function addComponentsToPath(testCase)
            components = fullfile(fileparts(mfilename("fullpath")), "..", "components");
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(components));
        end

        function hideFigures(testCase)
            original = get(groot, "DefaultFigureVisible");
            set(groot, "DefaultFigureVisible", "off");
            testCase.addTeardown(@set, groot, "DefaultFigureVisible", original);
        end
    end

    methods (Static)
        function wf = tone(power, bin, n, fs)
            % n samples at fs [Hz] of a tone of power [W] on DFT bin bin.
            wf = Waveform(sqrt(power) * exp(2j * pi * bin * (0:n - 1).' / n), fs, 24.125e9);
        end

        function wf = withSamples(wf, samples)
            % Assign samples as a component would.
            wf.Samples = samples;
        end
    end

    %% Signal
    methods (Test)
        function storesSignal(testCase)
            wf = Waveform([1 2j -3], 4e3, 24.125e9, Label="Test");
            testCase.verifyEqual(wf.Samples, [1; 2j; -3]);
            testCase.verifyEqual(wf.SampleRate, 4e3);
            testCase.verifyEqual(wf.CarrierFrequency, 24.125e9);
            testCase.verifyEqual(wf.Label, "Test");
            testCase.verifyEqual(wf.State, "conducted");
            testCase.verifyEqual(wf.Time, [0; 0.25e-3; 0.5e-3]);
            testCase.verifyEqual(wf.Power, 14 / 3, RelTol=1e-15);
            testCase.verifyEqual(wf.Wavelength, 299792458 / 24.125e9, RelTol=1e-15);
            unnamed = Waveform(1, 1, 0);
            testCase.verifyEqual(unnamed.Label, "Waveform");
        end

        function componentsCanReplaceSamples(testCase)
            wf = Waveform(ones(4, 1), 1e3, 0);
            wf.Samples = [1 2];
            testCase.verifyEqual(wf.Samples, [1; 2]);
            testCase.verifyEqual(wf.Time, [0; 1e-3]);
        end

        function allArgumentsAreRequired(testCase)
            args = {ones(4, 1), 1e3, 24e9};
            for n = 0:numel(args) - 1
                testCase.verifyError(@() Waveform(args{1:n}), "MATLAB:minrhs");
            end
        end

        function invalidSignalsThrow(testCase)
            testCase.verifyError(@() Waveform(ones(4), 1e3, 0), "MATLAB:validation:IncompatibleSize");
            testCase.verifyError(@() Waveform([], 1e3, 0), "MATLAB:validators:mustBeNonempty");
            testCase.verifyError(@() Waveform([1 NaN], 1e3, 0), "MATLAB:validators:mustBeFinite");
            testCase.verifyError(@() Waveform(1, 0, 0), "MATLAB:validators:mustBePositive");
            testCase.verifyError(@() Waveform(1, 1e3, -1), "MATLAB:validators:mustBeNonnegative");
            testCase.verifyError(@() WaveformTest.withSamples(Waveform(ones(4, 1), 1e3, 0), ones(4)), ...
                "MATLAB:validation:IncompatibleSize");
        end
    end

    %% Spectrum
    methods (Test)
        function toneFallsInItsBin(testCase, numSamples)
            fs = 8e3;
            wf = WaveformTest.tone(0.01, 5, numSamples, fs);
            [P, f] = wf.spectrum();
            testCase.verifySize(P, [numSamples 1]);
            testCase.verifyEqual(diff(f), fs / numSamples + zeros(numSamples - 1, 1), RelTol=1e-12);
            testCase.verifyEqual(f(floor(numSamples / 2) + 1), 0, "Zero frequency is in the middle.");
            [peak, i] = max(P);
            testCase.verifyEqual(f(i), 5 * fs / numSamples, RelTol=1e-12);
            testCase.verifyEqual(peak, 0.01, RelTol=1e-12);
        end

        function spectrumSumsToPower(testCase, numSamples)
            rng(1);
            wf = Waveform(complex(randn(numSamples, 1), randn(numSamples, 1)), 1e3, 0);
            testCase.verifyEqual(sum(wf.spectrum()), wf.Power, RelTol=1e-12);
        end
    end

    %% Plots
    methods (Test)
        function timePlotShowsInPhaseAndQuadrature(testCase)
            wf = WaveformTest.tone(1, 3, 200, 20e3);
            [ax, unit] = wf.plotTime();
            testCase.addTeardown(@close, ax.Parent);
            testCase.verifyEqual(unit, 1e-3);
            testCase.verifyEqual(string(ax.XLabel.String), "Time [ms]");
            i = findobj(ax, "Type", "line", "DisplayName", "I");
            q = findobj(ax, "Type", "line", "DisplayName", "Q");
            testCase.verifyEqual(i.XData(:), wf.Time / 1e-3);
            testCase.verifyEqual(i.YData(:), real(wf.Samples));
            testCase.verifyEqual(q.YData(:), imag(wf.Samples));
            testCase.verifyEqual(string(ax.Title.String), wf.Label);
        end

        function timePlotZoomsAndShowsRealSignalOnce(testCase)
            wf = Waveform(cos(2 * pi * (0:999).' / 50), 1e6, 0);
            ax = wf.plotTime(TimeLimits=[0 20e-6]);
            testCase.addTeardown(@close, ax.Parent);
            testCase.verifyNumElements(findobj(ax, "Type", "line"), 1);
            testCase.verifyEqual(string(ax.XLabel.String), "Time [\mus]");
            testCase.verifyEqual(ax.XLim, [0 20], AbsTol=1e-12);
        end

        function spectrumPlotShowsPowerAroundCarrier(testCase)
            wf = WaveformTest.tone(0.01, 5, 64, 8e3);
            [ax, unit] = wf.plotSpectrum(DynamicRange=60);
            testCase.addTeardown(@close, ax.Parent);
            [P, f] = wf.spectrum();
            line = findobj(ax, "Type", "line");
            testCase.verifyEqual(unit, 1e3);
            testCase.verifyEqual(line.XData(:), f / 1e3);
            testCase.verifyEqual(line.YData(:), max(10 * log10(P / 1e-3), -50));
            testCase.verifyEqual(ax.YLim, [-50 15], AbsTol=1e-9);
            testCase.verifyEqual(ax.XLim, [-4 4]);
            testCase.verifyEqual(string(ax.XLabel.String), "Frequency - 24.125 GHz [kHz]");
            testCase.verifySubstring(ax.Subtitle.String, "10.0 dBm");

            baseband = Waveform(ones(64, 1), 8e3, 0);
            ax = baseband.plotSpectrum();
            testCase.addTeardown(@close, ax.Parent);
            testCase.verifyEqual(string(ax.XLabel.String), "Frequency [kHz]");
            line = findobj(ax, "Type", "line");
            testCase.verifyEqual([min(line.YData) max(line.YData)], [-50 30], ...
                "AbsTol", 1e-9, "Empty bins around a lone tone sit on the floor of the plot.");
        end

        function plotsIntoGivenAxes(testCase)
            fig = figure;
            testCase.addTeardown(@close, fig);
            layout = tiledlayout(fig, 1, 2);
            left = nexttile(layout);
            right = nexttile(layout);
            wf = WaveformTest.tone(1, 3, 64, 8e3);
            testCase.verifySameHandle(wf.plotTime(Parent=left), left);
            testCase.verifySameHandle(wf.plotSpectrum(Parent=right), right);
            testCase.verifyNumElements(findobj(fig, "Type", "axes"), 2);
        end

        function plotShowsEveryStage(testCase)
            stages = [WaveformTest.tone(1, 3, 64, 8e3), Waveform(ones(32, 1), 8e3, 0, Label="Second")];
            fig = plot(stages, TimeLimits=[0 2e-3], DynamicRange=40);
            testCase.addTeardown(@close, fig);
            axs = findobj(fig, "Type", "axes");
            testCase.verifyNumElements(axs, 4);
            titles = arrayfun(@(ax) string(ax.Title.String), axs);
            testCase.verifyEqual(sort(titles), ["Second"; "Second"; "Waveform"; "Waveform"]);
            testCase.verifyEqual(fig.Name, 'Waveform, Second');
        end
    end
end
