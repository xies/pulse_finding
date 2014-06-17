
[f,nC] = estimate_simulation_params(fits_wt,cells_wt);

%% Perform permultation analysis

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

%% Save the XYT coordinates of the randomized pulses

[phat,pci] = gamfit([freqOI{:}]);
frequency.fun = @(x) gamcdf(x,phat(1),phat(2));

pc_emp = [cellsOI.get_curated.num_fits];
N_pulse_count = hist(pc_emp,0:20);
N_pulse_count = cumsum(N_pulse_count)/sum(N_pulse_count);

p = cell(1,Nboot);
seq = cell(1,Nboot);

for i = 1:Nboot
    
    [p{i},seq{i}] = simulate_pulsing( ...
        cellsOI,fitsOI,frequency,N_pulse_count);
    
    display(['Done with ' num2str(i)]);
    
    for embryoID = 1:5
        
        this_p = p{i}([p{i}.embryoID] == embryoID);
        
        if embryoID == 1, embryoID = 8; end
        
        fIDs = cat(1,this_p.fitID);
        cy = cat(1,this_p.centroid_y);
        cx = cat(1,this_p.centroid_x);
        ct = cat(1,this_p.center);
        l = cat(1,this_p.cluster_label);
        
        M = cat(2,fIDs,cx,cy,ct,l);
        path = ['~/Desktop/Pulse xyt csv/Embryo ' ...
            num2str(embryoID) '/simulated/emb' num2str(embryoID) '_N' num2str(i) '.csv'];
        
%         csvwrite(path,M);
        
    end
end
