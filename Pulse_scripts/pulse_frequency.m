%% Pulse frequency
% Wildtype

fitsOI = fits.get_embryoID([11 12]);

[freqOI,centerOI] = cells.get_frequency(fitsOI);

%% WT plots

bins = linspace(0,300,30);
[freq_wt,center_wt] = cells.get_frequency(fits_wt);
% figure(1)
% [N_wt,bins] = hist( [freq_wt{:}], bins);
% bar( bins, N_wt/sum(N_wt) );
% xlim([0 300])
% xlabel('Time between pulses (sec)')
% ylabel('Probability')
% title('Wild-type')

% figure(1)
% subplot(2,1,1)
scatter([center_wt{:}], [freq_wt{:}], 100, 'filled', 'r')
xlim([-300 800])
xlabel('Developmental time (sec)');
ylabel('Time between pulses (sec)');
title('twist-RNAi injection')

%% twist

fits_incell = cellfun(@fits.get_fitID,...
    {cells.get_embryoID(6:10).fitID},'UniformOutput',0);

fits_label_incell = cell(1,numel(fits_incell));
fits_center_incell = cell(1,numel(fits_incell));
centers = cell(1,numel(fits_incell));
nnear = cell(1,numel(fits_incell));

for i = 1:numel(fits_incell)
    
    fits_incell{i} = fits_incell{i}.sort('center');
    fits_center_incell{i} = [fits_incell{i}.center];
    foo = [fits_incell{i}.cluster_label];
    fits_label_incell{i} = foo(1:end-1);
    foo = [fits_incell{i}.center];
    centers{i} = foo(1:end-1);
    foo = cat(1,fits_incell{i}.nearIDs);
    if ~isempty(foo)
        foo = cellfun(@numel,foo(:,6));
        nnear{i} = foo(1:end-1)';
    end
    
end

freq_twist = cellfun(@(x) diff(sort(x)), fits_center_incell, 'UniformOutput',0);
labels_twist = fits_label_incell;
center_twist = cellfun(@sort_pair_mean, fits_center_incell, 'UniformOutput',0);
[N_twist,bins] = hist([freq_twist{:}],bins);

%%

fits_incell = cellfun(@fits.get_fitID,...
    {cells.get_embryoID(14).fitID},'UniformOutput',0);

fits_label_incell = cell(1,numel(fits_incell));
fits_center_incell = cell(1,numel(fits_incell));
centers = cell(1,numel(fits_incell));
nnear = cell(1,numel(fits_incell));

for i = 1:numel(fits_incell)
    
    fits_incell{i} = fits_incell{i}.sort('center');
    fits_center_incell{i} = [fits_incell{i}.center];
    foo = [fits_incell{i}.cluster_label];
    fits_label_incell{i} = foo(1:end-1);
    foo = [fits_incell{i}.center];
    centers{i} = foo(1:end-1);
    foo = cat(1,fits_incell{i}.nearIDs);
    if ~isempty(foo)
        foo = cellfun(@numel,foo(:,6));
        nnear{i} = foo(1:end-1)';
    end
end

freq_cta = cellfun(@(x) diff(sort(x)), fits_center_incell, 'UniformOutput',0);
center_cta = cellfun(@sort_pair_mean, fits_center_incell, 'UniformOutput',0);

scatter([center_cta{:}],[freq_cta{:}]);
xlabel('Developmental time (sec)');

%% cta (seperate two cta populations)... see cta_clustering.m

cta_cells8 = cells( [cells.embryoID] == 8 );
cta_cells9 = cells( [cells.embryoID] == 9 );

fits_incell8 = cellfun(@fits_all.get_fitID, ...
    {cells(  c8([cta_cells8.cellID]) == 1 ).fitID}, 'UniformOutput', 0);
fits_incell9 = cellfun(@fits_all.get_fitID, ...
    {cells(  c9([cta_cells9.cellID]) == 1 ).fitID}, 'UniformOutput', 0);
fits_incell = [fits_incell8, fits_incell9];

fits_center_incell = cell(1,numel(fits_incell));
for i = 1:numel(fits_incell)
    fits_center_incell{i} = [fits_incell{i}.center];
end

freq_cta1 = cellfun(@(x) diff(sort(x)), fits_center_incell, 'UniformOutput',0);

fits_incell8 = cellfun(@fits_all.get_fitID, ...
    {cells(  c8([cta_cells8.cellID]) == 2 ).fitID}, 'UniformOutput', 0);
fits_incell9 = cellfun(@fits_all.get_fitID, ...
    {cells(  c9([cta_cells9.cellID]) == 2 ).fitID}, 'UniformOutput', 0);
fits_incell = [fits_incell8, fits_incell9];

fits_center_incell = cell(1,numel(fits_incell));
for i = 1:numel(fits_incell)
    fits_center_incell{i} = [fits_incell{i}.center];
end

freq_cta2 = cellfun(@(x) diff(sort(x)), fits_center_incell, 'UniformOutput',0);

%% Plot results as histograms

bins = linspace(0,250,50);

[N_wt, bins] = hist( [freq_wt{:}] , bins);
[N_twist,bins] = hist( [freq_twist{:}], bins);
% [N_cta1, bins] = hist( [freq_cta1{:}] , bins);
% [N_cta2, bins] = hist( [freq_cta2{:}] , bins);
% N_cta = N_cta1 + N_cta2;

bar( bins, ...
    cat( 1,N_wt/sum(N_wt), ...
    N_twist/sum(N_twist) ...
    )', 'Grouped');

set(gca,'XLim',[0 250]);

xlabel('Period between pulses (sec)')
ylabel('Probability')
legend(['Wild-type, N = ' num2str(sum(N_wt))], ...
    ['twist, N = ' num2str(sum(N_twist))] ...
    );

% figure, bar( bins, cat(1,N_cta1/sum(N_cta1), N_cta2/sum(N_cta2))' ,'Grouped');
% set(gca,'Xlim',[0 250]);
% xlabel('Period between pulses (sec)')
% ylabel('Probability')
% legend(['cta (constricting), N = ' num2str(sum(N_cta1))], ...
%     ['cta (expanding), N = ' num2str(sum(N_cta2))] ...
%     );

%% Count number of pulses per cell for different genotypes

N_wt = hist( [cells( [cells.embryoID] < 6 ).num_fits] , 1:15);
N_twist = hist( [cells( ismember([cells.embryoID], [6 7]) ).num_fits] ,1:15);
% N_cta = hist( [cells( [cells.embryoID] > 7 ).num_fits] ,1:15);

bar(1:15, ...
    cat( 1, N_wt, N_twist)' );

legend(['Wild-type, cells = ' num2str(sum(N_wt))], ...
    ['twist, cells = ' num2str(sum(N_twist))] ...
    );

%% Estimate frequency via MTrackJ data

mdf_mat = read_mdf('~/Desktop/SqhmCh SpiGFP alphaCat 5/alpha_cat_tracks');
I = [1 1 1 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 ...
    7 7 7 7 7 8 8 8 9 9 9 10 10 10 11 11 11 12 12 12 12 13 13 13];

alphacat_tracks = accumarray(mdf_mat(:,1),mdf_mat(:,end),[],@(x) {x});
center_cat = nonans(cellfun(@mean,alphacat_tracks));
centers_cat = accumarray(I',center_cat,[],@(x) {sort(x,'ascend')});
freq_cat = cellfun(@(x) diff(x)'*12.3,centers_cat,'UniformOutput',0);
