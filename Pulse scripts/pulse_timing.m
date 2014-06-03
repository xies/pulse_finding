% Pulse timing

embryoID = 1:5;
% embryoID = 6:10;
% embryoID = 11:13;

fitsOI = fits_bs.get_embryoID(embryoID);
x_limits = [-300 800];
bins = linspace(x_limits(1),x_limits(2),50);

%% by bin

colors = pmkmp(10);

clear N
for i = 1:10
    hold on;
%     Nwt(i,:) = hist([fitsOI([fitsOI.bin] == i).center],bins);
%     N(i,:) = hist([fitsOI([fitsOI.bin] == i).center],bins);
%     plot(bins,cumsum(N(i,:))/sum(N(i,:)),'Color',colors(i,:));
%     subplot(10,1,i)
%     fitsOI = fitsOI.bin_fits;
    N(i,:) = plot_pdf([fitsOI([fitsOI.bin] == i).center],bins,'FaceColor',colors(i,:));
%     xlim(x_limits);
    mean_bin_center(i,embryoID) = mean([fitsOI([fitsOI.bin]== i).center]);
    width(i,embryoID) = std([fitsOI([fitsOI.bin]== i).center]);
end
hold off

% Nwt = N;
Ntwist = N;

% imagesc(bins,1:10,N); colormap hot; axis xy;
ylabel('Probability')
xlabel('Developmental time (sec)');
% end

%% Heatmap instead of PDF/CDF line plots?

subplot(2,1,1);
imagesc(bins,1:10,Nwt);
colormap hot; colorbar;
axis xy

subplot(2,1,2);
imagesc(bins,1:10,Ntwist);
colormap hot; colorbar;
axis xy

%% KL div between one bin and the next?

bins = linspace(x_limits(1),x_limits(2),1000);
for i = 1:9
    this = [fitsOI([fitsOI.bin] == i).center];
    next = [fitsOI([fitsOI.bin] == i + 1).center];
    
    pThis = ksdensity( this, bins);
    pNext = ksdensity( next, bins);
    
    KL(i) = kldiv(bins,pThis/sum(pThis),pNext/sum(pNext));
    
end
plot(KL)
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
N = zeros(2,numel(bins));
for i = 1:num_clusters
    N(i,:) = hist([fitsOI([fitsOI.cluster_label]==i).center],bins);
end

% N = bsxfun(@rdivide, N, sum(N,2));

h = plot(bins,bsxfun(@rdivide,N,sum(N,2))');
% for i = 1:num_clusters
%     set(h(i),'FaceColor',colors{i});
% end
legend(behaviors{:})
xlabel('Developmental time (sec)')
ylabel('Probability')
xlim(x_limits)

%% behavior by temporal bins

left = [-Inf	0];
right = [0  Inf];
N = zeros( numel(left), 6);

for i = 1:numel(left)
    
    filter = @(x) ([x.center] > left(i) & [x.center] <= right(i));
    N(i,:) = hist( [fitsOI(filter(fitsOI)).cluster_label], 1:6);
    
end
%filter out non-clustered pulses
N(:,num_clusters+1) = [];

h = bar(1:numel(left), bsxfun(@rdivide,N,sum(N,2)), 'stacked','LineStyle','none');

for i = 1:2
    set(h(i),'FaceColor',colors{i});
end
legend([behaviors,'N/A']);
ylabel('Count');
