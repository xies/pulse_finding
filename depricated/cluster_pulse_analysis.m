%% Stability analysis of FCM clustering

[RI,randRI] = pulse.fcm_stabiliy(Ks2try);

%%

num_clusters = 3;

fits.fcm_cluster(num_clusters,'corrected_area_norm',3);

%%

clear cluster*

for i = 1:num_clusters
    
    eval(['cluster' num2str(i) ' = fits([fits.cluster_label] == ' num2str(i) ');']);
    
end

behaviors = {'Ratcheted',...
    'Un-ratcheted', ...
    'Unconstricting'};

colors = {'b','m','r'};

%%

figure

% fits_cta = fits.get_embryoID( 8 );
% [N_const,bins] = hist(revorder([fits_cta( c8([fits_cta.cellID]) == 1 ).cluster_label]), ...
%     1:num_clusters);
% [N_exp,bins] = hist(revorder([fits_cta( c8([fits_cta.cellID]) == 2 ).cluster_label]), ...
%     1:num_clusters);

[N_wt] = hist( [fits_wt.cluster_label], 1:num_clusters);
[N_twist] = hist( [fits_twist.cluster_label], 1:num_clusters);
[N_cta] = hist( [fits_cta.cluster_label], 1:num_clusters);

N_wt = N_wt(order);
N_twist = N_twist(order);
N_cta = N_cta(order);

bar(1:3, ...
    cat(1,N_wt/sum(N_wt),N_twist/sum(N_twist),N_cta/sum(N_cta)),'stacked')
set(gca,'XTickLabel',{'Wild-type','twist','cta'});
% bar(1:2,cat(1,N_const/sum(N_const),N_exp/sum(N_exp)),'stacked')
% set(gca,'XTickLabel',{'Constricting','Expanding'})
ylabel('Probability')
legend(entries{:});

%% summary of clusters

figure
for i = 1:num_clusters
    
    eval(['this_cluster = cluster' num2str(i) '_wt.sort(''cluster_weight'');']);
    cluster_area = cat(1,this_cluster.corrected_area_norm);
    
    subplot(2,num_clusters,i);
    [X,Y] = meshgrid( fits(1).corrected_time,1:numel(this_cluster) );
    pcolor( X,Y, cluster_area );
    shading flat, caxis([-7 7]),colorbar;
        title(['Cluster ' num2str(i) ]);
    xlabel('Pulse time (sec)')
    
    subplot(2,num_clusters,i+num_clusters);
    weights = cat(1, this_cluster.cluster_weight);
    weights = max(weights,[],2);
    shadedErrorBar( fits(1).corrected_time, ...
        nanwmean(cluster_area,weights), nanstd(cluster_area) , colors{i});
    xlabel('Pulse time (sec)')
    ylim( [-6 6] );
    set(gca,'XTick',[-40 0 40]);
    
end

% figure
% % 2D histogram heatmap view
% for i = 1:num_clusters
%     subplot(1,num_clusters,i);
%     eval(['A = cat(1,cluster' num2str(i) '_wt.corrected_area_norm);']);
%     bins = -5:.1:5;
%     N = hist(A,bins);
%     imagesc(fits(1).corrected_time,bins,bsxfun(@rdivide,N,sum(N)));
%     axis square xy;
%     colormap hot;
%     colorbar;
%     xlabel('Pulse time (sec)')
% end

%% Breakdown behavior by embryoID

clear N

N(1,:) = hist([fits_wt.cluster_label],1:4);
N(2,:) = hist([fits_twist.cluster_label],1:4);
N(3,:) = hist([fits_control.cluster_label],1:4);

% N(:,4) = [];

bar(bsxfun(@rdivide, N, sum(N,2)) ,'stacked' );
xlim([0 4])
set(gca,'XTickLabel',{'WT','twist','cta'});
% xlabel('EmbryoID')
legend(behaviors{:},'N/A');


