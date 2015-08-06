function fits_new = bootstrap_cluster_label(pulse)
% Perform intra-embryo bootstrapping of cluster labels

%this is a test I am adding a comment

fits = [pulse.fits];

% now i deleted a line

x=10
labels = zeros(1,numel(fits));

for i = embryoIDs
    
    this = fits.get_embryoID(i);
    l = cat(1,fits.cluster_label);
    labels( [fits.embryoID] == i ) = ...
        l( randperm(numel(this)) );
    
end

fits_new = fits.copy; % Deep copy

for i = 1:numel(fits)
    
    fits_new(i).cluster_label = labels(i);
    fits(i).bootstrapped = 1;

end

end % bootstrap_cluster_label

 pint " cats! "