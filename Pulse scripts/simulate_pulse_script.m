
[f,nC] = estimate_simulation_params(fits_wt,cells_wt);

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

%% Look at multiple permutation analysis

dir = '~/Desktop/';
filebase = 'wild_type.mat.iter_';
fileend = '_permutation.mat';

MC = cell(1,2);
for i = 1:5
    
    filename = [dir,filebase, num2str(i), fileend];
    data = load(filename);
    
    MC{i} = data.MC;
    fbs{i} = data.fits;
    cbs{i} = data.cells;
    
    tau = 6; % neighborhood time window
    clear opt temporal_bins
    temporal_bins(1,:) = [-Inf];
    temporal_bins(2,:) = [Inf];
    
    opt.normalize = 'off';
    opt.breakdown = 'off';
    opt.xlim = [2 4];
    i
%     zscores_wt(i,:) = plot_mc_results(MC{i},tau,temporal_bins,opt);
    
end

%%
