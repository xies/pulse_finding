%%

dir = '~/Dropbox (MIT)/Adam/Pulse figures/MCstackID/wild_type/simulated_permute/';
filebase = 'wild_type.mat.iter_';
fileend = '_permutation.mat';

Nsim = 10;

%%

MC2plot = cell(1,Nsim);
for i = 1:Nsim
    
    path = [dir, filebase, num2str(i), fileend];
    data = load(path);
    MC2plot{i} = data.MC;
%     keyboard
    fbs{i} = data.fits;
    cbs{i} = data.cells;
    i
end

%% Look at multiple permutation analysis

for i = 1:10
    
    tau = 6; % neighborhood time window
    clear opt temporal_bins
    temporal_bins(1,:) = [-Inf];
    temporal_bins(2,:) = [Inf];
    
    opt.normalize = 'off';
    opt.breakdown = 'off';
    
    opt.xlim = [2 4];
    i
    zscores_sim(i,:) = plot_mc_results(MC2plot{i},tau,temporal_bins,opt);
    
end

%%

