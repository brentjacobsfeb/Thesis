classdef phasor
    properties
        amplitude = 1
        frequency = 1
        phase = 0
    end

    methods
        function obj = phasor(amplitude, frequency, phase)
            if nargin >= 1
                obj.amplitude = amplitude;
            end
            if nargin >= 2
                obj.frequency = frequency;
            end
            if nargin >= 3
                obj.phase = phase;
            end
        end
    end
end