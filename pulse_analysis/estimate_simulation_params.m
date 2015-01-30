function [freq,neighbor_count] = estimate_simulation_params(fits,cells)
%ESTIMATE_SIMULATION_PARAMETERS Returns frequency and neighbor cell count
% distributions for a given set of pulses (FITTED).
%
% USAGE: [freq,nCount] = estimate_simulation_params(fits,cells)


% Check inputs
if numel([fits.cluster_label]) ~= numel(fits)
    error('All pulses must have a cluster label');
end

% bins = linspace(0,300,30);
% 
% fits_incell = cellfun(@fits.get_fitID, {cells.fitID}, 'UniformOutput',0);
% fits_center_incell = cell(1,numel(fits_incell));
% 
% for i = 1:numel(fits_incell)
%     fits_incell{i} = fits_incell{i}.sort('center');
%     fits_center_incell{i} = [fits_incell{i}.center];
% end

f = cells.get_frequency(fits);

[phat,pci] = gamfit([f{:}]);
freq.fun = @(x) gamcdf(x,phat(1),phat(2));

% estimate nc
num_neighbors = zeros(1,numel(fits));
index = 0;

for i = unique([fits.embryoID]);
    
    f = fits.get_embryoID(i);
    
    N = cells.get_embryoID(i).get_adjacency_matrix;
    
    for j = 1:numel(f)
        
        index = index + 1;
        this_fit = f(j);
        this_conn = N(this_fit.cellID,:,this_fit.center_frame);
        
        num_neighbors(index) = numel(this_conn(this_conn > 0));
        
    end
    
end

neighbor_count = cell(1,max([fits.cluster_label]));
for i = 1:max([fits.cluster_label])
    neighbor_count{i} = ...
        hist(num_neighbors([fits.cluster_label] == i),0:10) ...
        /numel(find([fits.cluster_label]==i));
end

end