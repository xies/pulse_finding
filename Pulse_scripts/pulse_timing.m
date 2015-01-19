% Pulse timing

embryoID = 16;

fitsOI = fits.get_embryoID(embryoID);
x_limits = [-200 400];
bins = linspace(x_limits(1),x_limits(2),50);

% by bin

colors = pmkmp(10);

clear N
for i = 1:10
    hold on;
%     Nwt(i,:) = hist([fitsOI([fitsOI.bin] == i).center],bins);
    N(i,:) = hist([fitsOI([fitsOI.bin] == i).center],bins);
    plot(bins,cumsum(N(i,:))/sum(N(i,:)),'Color',colors(i,:));
    
    fitsOI = fitsOI.bin_fits;
%     N(i,:) = plot_pdf([fitsOI([fitsOI.bin] == i).center],bins,'FaceColor',colors(i,:));
    xlim(x_limits);
    mean_bin_center(i,embryoID) = mean([fitsOI([fitsOI.bin]== i).center]);
    width(i,embryoID) = std([fitsOI([fitsOI.bin]== i).center]);
end
hold off

Nwt = N;
% Ntwist = N;
% 
% imagesc(bins,1:10,N); colormap hot; axis xy;
% ylabel('Probability')
% xlabel('Developmental time (sec)');
% end

% Heatmap instead of PDF/CDF line plots?

% subplot(2,1,1);
% subplot(2,3,embryoID-5)
% 
% imagesc(bins,1:10,bsxfun(@rdivide,Nwt,sum(Nwt,2)) );
% title(['Embryo ' num2str(embryoID)])
% colormap hot; colorbar;
% axis xy


% subplot(2,1,2);
% imagesc(bins,1:10,bsxfun(@rdivide,Ntwist,sum(Ntwist,2)) );
% colormap hot; colorbar;
% axis xy

%% JS div between amplitude bins

bins = linspace(x_limits(1),x_limits(2),1000);

JSD_wt = zeros(10); JSD_twist = zeros(10);
for i = 1:10
    for j = 1:10
        
        % Use KDE to estimate
%         this = [fitsOI([fitsOI.bin] == i).center];
%         that = [fitsOI([fitsOI.bin] == j).center];
%         pThis = ksdensity( this, bins);
%         pThat = ksdensity( that, bins);        
%         JSD(i,j) = kldiv(bins,pThis/sum(pThis),pThat/sum(pThat),'js');

        % Use histogram estimation
        JSD_wt(i,j) = kldiv( linspace(x_limits(1),x_limits(2),size(Nwt,2)), ...
            (Nwt(i,:) + eps)/sum(Nwt(i,:)), ...
            (Nwt(j,:) + eps)/sum(Nwt(j,:)) , 'js');
        
        JSD_twist(i,j) = kldiv( linspace(x_limits(1),x_limits(2),size(Ntwist,2)), ...
            (Ntwist(i,:) + eps)/sum(Ntwist(i,:)), ...
            (Ntwist(j,:) + eps)/sum(Ntwist(j,:)) , 'js');
    end
end
% subplot(2,1,1);
% mesh(JSD_wt); colormap hot, caxis([0 max(JSD_wt(:))])
% subplot(2,1,2);
% mesh(JSD_twist); colormap hot, caxis([0 max(JSD_wt(:))])

[X,Y] = meshgrid(1:10);
surf( X,Y, triu(JSD_wt) + tril(JSD_twist) )

%% by behavior in CDF

colors = {'b-','m-','r-'};
for i = 1:num_clusters
    hold on
    plot_cdf([fitsOI.get_cluster(i).center],bins,colors{i});
    xlim(x_limits)
end

xlabel('Developmental time (sec)')
ylabel('Probability')
legend(behaviors{:})

%% behavior in count/PDF

colors = {'b','m','r'};
N = zeros(2,numel(bins)); means = zeros(1,num_clusters);
for i = 1:num_clusters
    x = [fitsOI([fitsOI.cluster_label]==i).center];
    N(i,:) = hist(x,bins);
    means(i,:) = nanmean(x);
end

% N = bsxfun(@rdivide, N, sum(N,2));

h = plot(bins,bsxfun(@rdivide,N,sum(N,2))');
for i = 1:num_clusters
    hold on,
    ylim([0 0.2])
    vline(means(i,:),colors{i});
    set(h(i),'Color',colors{i});
end
% legend(behaviors{:})
xlabel('Developmental time (sec)')
ylabel('Probability')
xlim(x_limits)

%% behavior by temporal bins

[~,which_bin] = histc([fitsOI.center],bins);

clear N
for i = 1:numel(bins)
    N(i,:) = hist( [fitsOI( which_bin == i ).cluster_label], 1:4);
end

%filter out non-clustered pulses
N(:,4) = [];
bar(bins,bsxfun(@rdivide,N(:,1:3),sum(N(:,1:3),2)),'stacked')
xlim(x_limits)



