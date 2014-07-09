single_bin = fits_wt([fits_wt.bin] == 9);

nearIDs = cat(1,single_bin.nearIDs);
num_near = cellfun(@(x) numel(x(~isnan(x))), nearIDs);


for i = 1:numel(single_bin)
    p = polyfit(fits(1).corrected_time, single_bin(i).corrected_area_norm, 1);
    coeffs(i,:) = p;
end

for i = 1:10
    foo{i} = coeffs(num_near(:,6) == i,1);
end

[foo{cellfun(@isempty,foo)}] = deal(NaN);

[X,G] = make_boxplot_args(foo{:});

subplot(5,1,1:4)
boxplot(X,G)
xlabel('# of neighboring pulses')
ylabel('Slope of local area change')

subplot(5,1,5)
bar(cellfun(@numel,foo))
