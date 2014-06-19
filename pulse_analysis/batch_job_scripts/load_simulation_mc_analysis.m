%%

for i = 1:10
    path = ['~/Desktop/simulated pulses/wild_type.mat.iter_' num2str(i) '_permutation.mat'];
    data = load(path);
    MC2plot{i} = data.MC;
end

%% Select correct timing

for i = 1:10
    
    % select dataset
    MC = MC2plot{i};
    
    % MC = filter_mc(MC,ismember([fits_twist.embryoID],[ 10]));
    
    tau = 6; % neighborhood time window
    clear opt temporal_bins
    temporal_bins(1,:) = [-Inf];
    temporal_bins(2,:) = [Inf];
    
    opt.normalize = 'off';
    opt.breakdown = 'off';
    opt.xlim = [2 4];
    % opt.normalize = [5.06 5.00 5.29 5.01];
    
    zscores_wt(i,:) = plot_mc_results(MC,tau,temporal_bins,opt);
    
end

%%

