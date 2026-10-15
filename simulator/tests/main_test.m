%MAIN_TEST Run the complete measurement chain once with the default set-up.
%   Plots the antenna pattern, its footprint on the road and every stage of
%   the waveform, and prints the ground speed estimate. result holds the
%   antenna, the road, every waveform stage and the estimate (see help main).
close all;
clear;
addpath(fullfile(fileparts(mfilename("fullpath")), "..", "components"));

result = main();
