
filtered = fits(all(~isnan(corrected_area_norm),2));
X = cat(1,filtered.get_embryoID(11:12).corrected_area_norm);

%%

[coeff,score,latent,tsquared,explained,mu] = pca(X'); % covariance
% [coeff,score,latent,tsquared,explained,mu] = pca(zscore(X,1,2)'); % correlation

%% Plot variance explained v. PC number

plotyy(1:numel(explained),explained,1:numel(explained),cumsum(explained))

%% Do subspace projection

c3 = coeff(:,1:3);
labels = kmeans(zscore(X,1,2),3);

[~,order] = sort(labels);
imagesc(X(order,:)),caxis([-8 8])

rand_index(labels,[filtered.cluster_label]')

%%

colors = {'b','m','r','c','g'};
figure

% l = [filtered.cluster_label];
[~,l] = histc([filtered.embryoID],[0 5 10 13 14]);

for i = 1:max(l)
    scatter3(coeff(l==i,1),coeff(l==i,2),coeff(l==i,3),50,colors{i},'filled');
    hold on
end

legend('1','2','3','4');
xlabel('Ratcheted'),ylabel('Unratcheted'),zlabel('Unconstricting')
