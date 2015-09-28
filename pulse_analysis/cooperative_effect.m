%% Cooperative

embryoID = 1:5; txt = 'wt_4';
% embryoID = 6:10; txt = 'twist_4';
% embryoID = 11:15; txt = 'control_4';

fitsOI = pulse(embryoID).getFits;
cellsOI = pulse(embryoID).getCells;

time_windows = 20;
clear neighbor_defition
neighbor_defition.temporal.def = @(time_diff,tau) (abs(time_diff) < tau);
neighbor_defition.temporal.windows = time_windows;
neighbor_defition.spatial.def = 'identity';

pulse.find_near_fits(neighbor_defition);
nearIDs = cat(1,fitsOI.nearIDs);
num_near = cellfun(@(x) numel(x(~isnan(x))), nearIDs);

%% Export CVS to R

Mcr = nanmax( -diff(cat(1,fitsOI.corrected_area_norm),1,2) , [],2 );
mean_cr = nanmean( -diff( cat(1,fitsOI.corrected_area_norm),1,2) , 2);

range = nanmax( cat(1,fitsOI.corrected_area_norm),[],2 ) - ...
    nanmin( cat(1,fitsOI.corrected_area_norm),[],2 );

% csvwrite(['~/Dropbox (MIT)/' txt '.txt'], ...
%     cat(2,num_near(:,1),Mcr,range,[fitsOI.cluster_label]',[fitsOI.bin]') );

%%

clear medians
% bin = 10
for bin = 1:10
    
    single_bin = fitsOI([fitsOI.bin] == bin);
%     single_bin = single_bin([single_bin.cluster_label] == 1);
    
    if isempty(single_bin)
        continue
    end
    
    nearIDs = cat(1,single_bin.nearIDs);
    
    num_near = cellfun(@(x) numel(x(~isnan(x))), nearIDs);
    
    clear coeffs
    for i = 1:numel(single_bin)
        %         p = polyfit(fits(1).corrected_time, single_bin(i).corrected_area_norm, 1);
        p = nanmax(-diff(single_bin(i).corrected_area_norm));
        %         p = nanmean(-diff(Ca(i,:)));
        %         p = nanmax(single_bin(i).corrected_area_norm) - nanmin(single_bin(i).corrected_area_norm);
        coeffs(i) = p(1);
    end
    
    foo = cell(1,16);
    for i = 0:15
        foo{i+1} = coeffs(num_near(:,1) == i);
    end
    
    [foo{cellfun(@isempty,foo)}] = deal(NaN);
    % throw out only 1-element parts
    [foo{cellfun(@(x) numel(x) < 2,foo)}] = deal(NaN);
    
    figure(1)
    subplot(5,2,bin);
    boxplot(coeffs,num_near(:,1));
    hold on
    medians(:,bin) = cellfun(@nanmedian,foo);
    plot(medians(:,bin));
    ylim([.5 6])
    
end

% medians_wt = medians;
% medians_wt1 = medians;
% medians_wt2 = medians;
% medians_wt3 = medians;
% medians_twist = medians;
% medians_twist1 = medians;
% medians_twist2 = medians;
% medians_twist3 = medians;

%%

figure(2)
subplot(1,2,1),hold on
set(gca,'ColorOrder',pmkmp(10))
plot(0:15,medians_wt),xlim([0 10]),title('WT')
% ylim([-10 15])

subplot(1,2,2),hold all
set(gca,'ColorOrder',pmkmp(10))
plot(0:15,medians_twist),xlim([0 10]),title('twist')
% ylim([-10 15])

xlabel('# neighoring pulses (<60)'),ylabel('Mean centroid displacement')

%%

figure(1)
subplot(3,1,1),hold on
set(gca,'ColorOrder',pmkmp(10))
plot(0:15,medians_wt1),xlim([0 10]),title('WT')

subplot(3,1,2),hold all
set(gca,'ColorOrder',pmkmp(10))
plot(0:15,medians_wt2),xlim([0 10]),title('twist')

subplot(3,1,3),hold all
set(gca,'ColorOrder',pmkmp(10))
plot(0:15,medians_wt3),xlim([0 10]),title('twist')

xlabel('# neighoring pulses (<60)'),ylabel('Max constriction rate')

%%

figure(3)
subplot(3,1,1),hold on
set(gca,'ColorOrder',pmkmp(10))
plot(0:15,medians_twist1),xlim([0 10]),title('WT')

subplot(3,1,2),hold all
set(gca,'ColorOrder',pmkmp(10))
plot(0:15,medians_twist2),xlim([0 10]),title('twist')

subplot(3,1,3),hold all
set(gca,'ColorOrder',pmkmp(10))
plot(0:15,medians_twist3),xlim([0 10]),title('twist')

xlabel('# neighoring pulses (<60)'),ylabel('Max constriction rate')

%%

[X,G] = make_boxplot_args(foo{:});

subplot(5,1,1:4)
boxplot(X,G)
xlabel('# of neighboring pulses')
ylabel('Slope of local area change')

subplot(5,1,5)
bar(cellfun(@(x) numel(x(~isnan(x))),foo))
title('WT')
