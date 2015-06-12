function fits = bootstrap_cluster_label(fits)
% Perform intra-embryo bootstrapping of cluster labels
embryoIDs = unique([fits.embryoID]);
labels = zeros(1,numel(fits));

for i = embryoIDs
    
    this = fits.get_embryoID(i);
    l = cat(1,fits.cluster_label);
    labels( [fits.embryoID] == i ) = ...
        l( randperm(numel(this)) );
    
end

for i = 1:numel(fits)
    fits(i).cluster_label = labels(i);
    fits(i).bootstrapped = 1;
end

end % bootstrap_cluster_label