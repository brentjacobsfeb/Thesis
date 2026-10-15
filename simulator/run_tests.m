function results = run_tests(name)
%RUN_TESTS Run the simulator unit tests.
%   run_tests() runs every test in simulator/tests.
%   run_tests("RadarAntennaTest") runs one test file;
%   run_tests("RadarAntennaTest/taperLowersSidelobes") runs one test.
arguments
    name (1,1) string = ""
end

testsFolder = fullfile(fileparts(mfilename("fullpath")), "tests");
if name == ""
    suite = testsuite(testsFolder);
else
    suite = testsuite(fullfile(testsFolder, extractBefore(name + "/", "/") + ".m"));
    if contains(name, "/")
        suite = suite.selectIf("Name", "*" + extractAfter(name, "/") + "*");
    end
end

results = run(suite);
disp(table(results));
end
