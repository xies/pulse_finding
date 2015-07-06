function plot_aligned(cells,name2plot)

if nargin < 2, name2plot = 'area_sm'; end

embryoIDs = unique([cells.embryoID]);
color = hsv(numel(embryoIDs));

for i = 1:numel(embryoIDs)
    
    eID = embryoIDs(i);
    
    cells_in_embryo = cells.get_embryoID(eID);
    time = cells_in_embryo(1).dev_time;
    dt = mean(diff(time));
    data = cat(2,cells_in_embryo.(name2plot));
    
    H(i) = shadedErrorBar(time,...
        nanmean(data,2), nanstd(data,[],2), ...
        {'color',color(i,:)},1);

    labels{i} = ['Embryo ' num2str(eID) ', ' num2str(dt) ' sec/frame'];
    
    hold on
    
end

hold off
ylabel('Number of cells connected by myosin')
xlabel('Developmental time')
legend([H.mainLine],labels)

end