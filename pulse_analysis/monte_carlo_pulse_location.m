function MC = monte_carlo_pulse_location(pulse,opt)
%MONTE_CARLO_PULSE_LOCATION Generate spatially random pulse patterns by
% permutation analysis
% 
% USAGE: MC = monte_carlo_pulse_location(fits,cells,cx,cy,neighborID,opt)
%
% INPUT: fits - pulses
% 		 cells - cells
%		 neighborID - neighborhood of each cell
%		 opt - options for MC analysis
%		 	.Nboot - number of iterations
%			.time_windows - windows in time to check
%			.neighbor_def - neighborhood defintion
%			.savepath
%
% OUTPUT: MC
%
% xies@mit.edu Oct 2013

% Parse input options structure
Nboot = opt.Nboot;
time_windows = opt.timewindows;
neighbor_def = opt.neighbor_def;

% ---- Gather data for empirical pattern ----
pulse.find_near_fits(neighbor_def);
if strcmpi(opt.filter,'on')
    fitsOI = pulse.find_non_edge;
else
    fitsOI = pulse.fits;
end

nearIDs = cat(1,fitsOI.nearIDs);
% calculate num-neighbors for each fitsOI
num_near = cellfun(@(x) numel(x(~isnan(x))), nearIDs);
num_cells = cat(1,fitsOI.neighbor_cells);

% origin labels
labels = [fitsOI.cluster_label]';
% get all target labels
%         target_labels = ...
%             cellfun( @(x) [fits.get_fitID(x).cluster_label], ...
%             this_nearIDs,'UniformOutput',0);

% get centers
centers = [fitsOI.center];
CT = [ fitsOI(~isnan([fitsOI.nearest_neighbor])).center ];
nn = fitsOI.get_fitID([fitsOI.nearest_neighbor]);
dt = [nn.center] - CT;

% Should only record fits from non-edge of tissue cells
empirical.num_near = num_near;
empirical. angle
empirical.origin_labels = labels;
%         empirical.target_labels = target_labels;
empirical.centers = centers;
empirical.neighbor_windows = dt;
empirical.neighbor_cells = num_cells;
empirical.embryoID = [fitsOI.embryoID];

% ---- Gather data for simulated patterns ----

% Preallocate variables for simulated patterns
random_cell(Nboot).num_near = [];
random_cell(Nboot).origin_labels = [];
random_cell(Nboot).target_labels = [];
random_cell(Nboot).embryoID = [];
random_cell(Nboot).centers = [];
random_cell(Nboot).fitID = [];

for j = 1:Nboot
    
    tic
    % Simulate fit locations
    if strcmpi(opt.monte_carlo,'simulation')
        [f,~] = estimate_pulsing_params(pulse);
        pulse_sim = pulse.simulate_pulsing(f);
    else
        error('Can only use ''simulation'' optionl; others unstable')
    end
    % Find and update nearby fits
	pulse_sim.find_near_fits(neighbor_def);
    
    % Gather simulated fits
    if strcmpi(opt.filter,'on')
        fits_sim = pulse_sim.find_non_edge;
    else
        fits_sim = pulse_sim.fits;
    end
    
    nearIDs = cat(1,fits_sim.nearIDs);
    num_cells = cat(1,fits_sim.neighbor_cells);
    % tabulate num-neighbors for each fit
    num_near = cellfun(@(x) numel(x(~isnan(x))), nearIDs);
    % get origin labels
    labels = [fits_sim.cluster_label]';
    % break down all target labels (doesn't work anymore)
%     target_labels = ...
%         cellfun(@(x) [fits_bs.get_fitID(x).cluster_label], ...
%         nearIDs_cell,'UniformOutput',0);
    centers = [fits_sim.center];
    
    CT = [fits_sim(~isnan([fits_sim.nearest_neighbor])).center];
    nn = fits_sim.get_fitID([fits_sim.nearest_neighbor]);
    
    dt = [nn.center] - CT;
    
    random_cell(j).num_near = num_near;
    random_cell(j).origin_labels = labels;
%     random_cell(j).target_labels = target_labels;
    random_cell(j).centers = centers;
    random_cell(j).fitID = [fits_sim.fitID];
    random_cell(j).neighbor_windows = dt;
    random_cell(j).neighbor_cells = num_cells;
    random_cell(j).embryoID = [fits_sim.embryoID];
    
    T = toc;
    display(['Done with ' num2str(j) ' in ' num2str(T) ' seconds.']);
    
end

MC.empirical = empirical;
MC.random_cell = random_cell;
MC.neighbor_def = neighbor_def;
MC.time_windows = time_windows;
MC.option = opt;

if isfield(opt,'savepath')
    if ~isempty(opt.savepath)
        save(opt.savepath,'MC','fits','cells');
        display(['Saved to: ' opt.savepath])
    end
end

end
