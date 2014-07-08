%BATCH_RANDOMIZE
% To be used from a master shell script for use in batch jobs.

function batch_simulate_pulsing(INPUT_MAT_FILE,OUT_FILENAME)

[DIR_STR,name,ext] = fileparts(INPUT_MAT_FILE);

if ~strcmpi(ext,'.mat'),
	error('Needs a .mat file as input.');
end

disp(['Loading input pulse file: ' INPUT_MAT_FILE]);
load(INPUT_MAT_FILE);

time_windows = 10:10:100; % seconds

neighbor_defition.temporal.def = @(time_diff,tau) (time_diff < tau & time_diff > 0);
neighbor_defition.temporal.windows = time_windows;
neighbor_defition.spatial.def = 'identity';

fits = fits.find_near_fits(cells,neighbor_defition);
o.Nboot = 100;
o.timewindows = time_windows;
o.neighbor_def = neighbor_defition;
o.monte_carlo = 'simulation';
o.filter = 'off';

o.savepath = [OUT_FILENAME]

disp('Randomization analysis...')
MC = monte_carlo_pulse_location(fits,cells, o);

end
