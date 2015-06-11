
[f,nC] = estimate_simulation_params(fits,cells);

fitsOI = fits.get_embryoID(1:2);
cellsOI = cells.get_embryoID(1:2);

%% Save the XYT coordinates of the randomized pulses

[phat,pci] = gamfit([freq_wt{:}]);
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
    
    foo = cellfun(@diff, seq{i}, 'UniformOutput', 0);
    freq_sim{i} = [foo{:}];
    
%     for embryoID = 1:5
%         
%         this_p = p{i}([p{i}.embryoID] == embryoID);
%         
%         if embryoID == 1, embryoID = 8; end
%         
%         fIDs = cat(1,this_p.fitID);
%         cy = cat(1,this_p.centroid_y);
%         cx = cat(1,this_p.centroid_x);
%         ct = cat(1,this_p.center);
%         l = cat(1,this_p.cluster_label);
        
%         M = cat(2,fIDs,cx,cy,ct,l);
%         path = ['~/Desktop/Pulse xyt csv/Embryo ' ...
%             num2str(embryoID) '/simulated/emb' num2str(embryoID) '_N' num2str(i) '.csv'];
        
%         csvwrite(path,M);
        
%     end
end
