%BATCH_SIMULATE_PULSING
% To be used from a master shell script for use in batch jobs.

function batch_simulate_pulsing(INPUT_MAT_FILE,OUT_FILENAME)

[DIR_STR,name,ext] = fileparts(INPUT_MAT_FILE);

if ~strcmpi(ext,'.mat'),
	error('Needs a .mat file as input.');
end

disp(['Loading input pulse file: ' INPUT_MAT_FILE]);
load(INPUT_MAT_FILE);

% Estimate pulsing parameters
[f,nC] = estimate_simulation_params(fits,cells);

disp('Performing pulsing simulation...');
tic
[fits_bs,cells_bs] = ...
	fits.simulate_pulsing(cells,f);
T = toc;
disp(['Simulation finished in ' num2str(T) ' seconds'])

time_windows = 10:10:100; % seconds

neighbor_defition.temporal.def = @(time_diff,tau) (time_diff < tau & time_diff > 0);
neighbor_defition.temporal.windows = time_windows;
neighbor_defition.spatial.def = 'identity';

fitsOI = fits_bs.find_near_fits(cells_bs,neighbor_defition);
o.Nboot = 100;
o.timewindows = time_windows;
o.neighbor_def = neighbor_defition;
o.monte_carlo = 'permute';
o.filter = 'on';

o.savepath = [OUT_FILENAME]

disp('Permutation analysis...')
tic;
MC = monte_carlo_pulse_location(fitsOI,cells_bs, o);
T = toc;
disp(['Permutation analysis finished in ' num2str(T) ' seconds'])

end
