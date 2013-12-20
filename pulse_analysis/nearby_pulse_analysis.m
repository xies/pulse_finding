%% Nearby pulse analysis

fitsOI = fits_twist;

%%

time_windows = 10:10:100; % seconds

neighbor_defition = @(central, neighbors, tau) ...
    abs( [neighbors.center] - central.center ) < tau ...
    & ~( neighbors == central ) ...
    & ([neighbors.center] - central.center) < 0;

fitsOI = fitsOI.find_near_fits(time_windows,neighborID, neighbor_defition);

%%

nearIDs = cat(1,fitsOI.nearIDs);

% Convert to number of pulses
num_near = cellfun(@(x) numel(x(~isnan(x))), nearIDs);

%% MC stackID

entries = {'Ratcheted (stereotyped)','Ratcheted (weak)','Ratcheted (delayed)','Un-ratcheted','Stretched'};

o.Nboot = 50;
o.timewindows = time_windows;
o.savepath = ['~/Desktop/mc_stackID_wt_pcenter_Nboot' num2str(o.Nboot) '_k' num2str(num_clusters)];
o.neighbor_def = neighbor_defition;
MC_wt_pcenter = monte_carlo_pulse_location(fitsOI,cells,neighborID, o);

%% Select correct timing

% select dataset
MC = MC;

window = 3; % neighborhood time window
clear temporal_bins
temporal_bins(1,:) = [-Inf, -Inf 0];
temporal_bins(2,:) = [Inf,   0  Inf];

opt.breakdown = 'off';

plot_mc_results(MC,window,temporal_bins,opt)

%% Cross-validation Variance testing

MC = MC_wt_pcenter;

Nboot = numel(MC.random_cell);

window = 6;
Ncv = 1000;
subsamples = [10 50 75 100 250]; % out of Nboots

zscores_cell = zeros(numel(subsamples),100,5);
zscores_pulse = zeros(numel(subsamples),100,5);

for k = 1:numel(subsamples)
    for j = 1:Ncv
        
        IDs = randi(Nboot, 1, subsamples(k));
        random_cell = MC.random_cell( IDs );
        random_pulse = MC.random_pulse( IDs );
        
        num_emp = MC.empirical.num_near;
       
        num_cell = cat(3,random_cell.num_near);
        num_pulse = cat(3,random_pulse.num_near);
        
        labels_emp = MC.empirical.origin_labels;
        labels_cell = random_cell(1).origin_labels;
        labels_pulse = random_pulse(1).origin_labels;
        
        for i = 1:num_clusters
            
            % breakdown by cluster behavior
            this_count_emp = num_emp( labels_emp == i, window);
            this_count_cell = squeeze( ...
                num_cell( labels_cell == i, window, :) );
            this_count_pulse = squeeze( ...
                num_pulse( labels_pulse == i, window, :) );
            
            mean_of_cell = mean(this_count_cell,1);
            mean_of_pulse = mean(this_count_pulse,1);
            
            zscores_cell(k,j,i) = ...
                (mean(this_count_emp) - mean(mean_of_cell)) ./ std(mean_of_cell);
            zscores_pulse(k,j,i) = ...
                (mean(this_count_emp) - mean(mean_of_pulse)) ./ std(mean_of_pulse);
            
        end
        
    end
    k
end

Z = cat(4,zscores_cell,zscores_pulse);

for i = 1:5
    subplot(5,1,i)
    errorbar(subsamples,mean(Z(:,:,i,1),2),std(Z(:,:,i,1),[],2),'r-'),hold on
    errorbar(subsamples,mean(Z(:,:,i,2),2),std(Z(:,:,i,2),[],2),'g-')
    title(entries{i});
end