%% Cooperative

embryoID = 1:5;
% embryoID = 6:10;

fitsOI = fits.get_embryoID(embryoID).get_cluster(3);
cellsOI = cells.get_embryoID(embryoID);

clear neighbor_defition
neighbor_defition.temporal.def = @(time_diff,tau) (abs(time_diff) < tau & time_diff < 0);
neighbor_defition.temporal.windows = time_windows;
neighbor_defition.spatial.def = 'identity';

fitsOI = fitsOI.find_near_fits(cellsOI,neighbor_defition);

%%

clear medians
clear coeffs
% bin = 10
for bin = 1:10
    % figure(2)
    
    single_bin = fitsOI([fitsOI.bin] == bin);
    
    nearIDs = cat(1,single_bin.nearIDs);
    num_near = cellfun(@(x) numel(x(~isnan(x))), nearIDs);
    
    for i = 1:numel(single_bin)
%         p = polyfit(fits(1).corrected_time, single_bin(i).corrected_area_norm, 1);
        p = nanmax(-diff(single_bin(i).corrected_area_norm));
        coeffs(i) = p(1);
    end
    
    foo = cell(1,16);
    for i = 0:15
        foo{i+1} = coeffs(num_near(:,6) == i);
    end
    
    [foo{cellfun(@isempty,foo)}] = deal(NaN);
    
    medians(:,bin) = cellfun(@nanmean,foo);
    bin_coeff(bin) = nanmean(coeffs);
    
end

% medians_wt = medians;
% medians_wt1 = medians;
% medians_wt2 = medians;
medians_wt3 = medians;
% medians_twist = medians;
% medians_twist1 = medians;
% medians_twist2 = medians;
% medians_twist3 = medians;

%%

subplot(1,2,1),hold on
set(gca,'ColorOrder',pmkmp(10))
plot(0:15,medians_wt),xlim([0 10]),title('WT')

subplot(1,2,2),hold all
set(gca,'ColorOrder',pmkmp(10))
plot(0:15,medians_twist),xlim([0 10]),title('twist')

xlabel('# neighoring pulses (<60)'),ylabel('Max constriction rate')

%%
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

[X,G] = make_boxplot_args(foo{:});

subplot(5,1,1:4)
boxplot(X,G)
xlabel('# of neighboring pulses')
ylabel('Slope of local area change')

subplot(5,1,5)
bar(cellfun(@(x) numel(x(~isnan(x))),foo))
title('WT')
