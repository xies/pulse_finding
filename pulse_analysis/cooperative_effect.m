%% Cooperative

for bin = 1:10
    % figure(2)
    
    single_bin = fits_wt([fits_wt.bin] == bin);
    
    nearIDs = cat(1,single_bin.nearIDs);
    num_near = cellfun(@(x) numel(x(~isnan(x))), nearIDs);
    
    for i = 1:numel(single_bin)
%         p = polyfit(fits(1).corrected_time, single_bin(i).corrected_area_norm, 1);
        p = nanmax(diff(single_bin(i).corrected_area_norm));
        coeffs(i,bin) = p;
    end
    
    foo = cell(1,10);
    for i = 1:10
        foo{i} = coeffs(num_near(:,6) == i,bin);
    end
    
    [foo{cellfun(@isempty,foo)}] = deal(NaN);
    
    medians(bin,:) = cellfun(@nanmedian,foo);
    
end

%%

[X,G] = make_boxplot_args(foo{:});

subplot(5,1,1:4)
boxplot(X,G)
xlabel('# of neighboring pulses')
ylabel('Slope of local area change')

subplot(5,1,5)
bar(cellfun(@numel,foo))
title('WT')
