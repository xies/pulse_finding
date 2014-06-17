%BATCH_SIMULATE_PULSING
% To be used from a master shell script for use in batch jobs.

[f,nC] = estimate_simulation_params(fits,cells);

% Perform permultation analysis

for n = 1:50

    % generate randomized cells
    name = 'wt';

    fitsOI = fits.get_embryoID(1:5);
    cellsOI = cells.get_embryoID(1:5);

    [fits_bs,cells_bs] = fitsOI.simulate_pulsing(cellsOI,f,nC);

    fitsOI = fits_bs.get_embryoID(1:5);
    cellsOI = cells_bs.get_embryoID(1:5);

    time_windows = 10:10:100; % seconds
    
    clear neighbor_definition
    neighbor_defition.temporal.def = @(time_diff,tau) (time_diff < tau & time_diff > 0);
    neighbor_defition.temporal.windows = time_windows;
    neighbor_defition.spatial.def = 'identity';

    fitsOI = fitsOI.find_near_fits(cellsOI,neighbor_defition);

    nearIDs = cat(1,fitsOI.nearIDs);
    near_angles = cat(1,fitsOI.near_angles);

    % Convert to number of pulses
    num_near = cellfun(@(x) numel(x(~isnan(x))), nearIDs);

    % Generates permutation cells
    entries = {'Ratcheted (stereotyped)','Ratcheted (weak)','Ratcheted (delayed)','Un-ratcheted','Stretched'};

    clear o
    o.Nboot = 50;
    o.timewindows = time_windows;
    o.neighbor_def = neighbor_defition;
    o.monte_carlo = 'permute';
    o.filter = 'on';

    % o.savepath = [];
    o.savepath = ...
        ['~/Desktop/simulated pulses/mc_stackID_', name, '_', ...
        'iter_', num2str(n), '_', neighbor_defition.spatial.def, ...
        '_Nboot', num2str(o.Nboot), '_', o.monte_carlo, '_neighborfilt_', o.filter, ...
        '_k' num2str(num_clusters)];

    MC_wt_sim{n} = monte_carlo_pulse_location(fitsOI,cellsOI, o);

end