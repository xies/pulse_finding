% Pulse strength

fitsOI = fits.get_embryoID(1:5);
fitsOI = fitsOI.bin_fits;

%%

N = zeros(10,num_clusters);
for i = 1:num_clusters
    
    N(:,i) = hist([fitsOI.get_cluster(i).bin]);
    
end

N = bsxfun(@rdivide,N,sum(N))';

figure(1)
bar(N,'stacked');
set(gca,'XTickLabel',behaviors);
ylabel('Probability')
ylim([0 1])
colormap(pmkmp(10))

figure(200)
plot(cumsum(N,2)');
legend(behaviors{:});
ylabel('Cumulative probability')

%% Local rank

nearIDs = cat(1,fitsOI.nearIDs);
nearIDs = nearIDs(:,2);
rank = zeros(1,numel(fitsOI));
for i = 1:numel(fitsOI)
    nearby_amps = [fitsOI.get_fitID(nearIDs{i}).amplitude];
    foo = tiedrank([fitsOI(i).amplitude nearby_amps]);
    rank(i) = foo(1);
end
