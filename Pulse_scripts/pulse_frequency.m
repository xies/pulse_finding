%% Pulse frequency
% Wildtype

fitsOI = fits.get_embryoID(1);

[freqOI,centerOI] = cells.get_frequency(fitsOI);
xlimits = [-300 800];
ylimits = [0 500];

%% Plots

bins = linspace(0,300,30);

% Histogram of freq distribution
figure(2)
subplot(2,1,1)
[N_wt,bins] = hist( [freqOI{:}], bins);
plot( bins, N_wt/sum(N_wt) );
hold on, vline(mean([freqOI{:}]))
xlim([0 300])
xlabel('Time between pulses (sec)')
ylabel('Probability')
title('char-RNAi')

% Scatter plot of dynamic freq change
figure(2)
subplot(2,1,2)
scatter([centerOI{:}], [freqOI{:}], 100, 'filled', 'r');

% get trendline
p = polyfit([centerOI{:}], [freqOI{:}],1);
hold on, xi = linspace(xlimits(1),xlimits(2),1000);
plot(xi,polyval(p,xi),'k-')

xlim(xlimits); ylim(ylimits);
xlabel('Developmental time (sec)');
ylabel('Time between pulses (sec)');
title('Wild-type')

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
