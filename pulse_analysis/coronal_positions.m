coronal_cx = cells.get_corona_measurement(centroids_x);
coronal_cy = cells.get_corona_measurement(centroids_y);

[T,C] = size(coronal_cx);

%%

ccx_centered = cellfun(@minus, coronal_cx, ...
    mat2cell(centroids_x, ones(1,T), ones(C,1)), ...
    'UniformOutput', 0);
ccy_centered = cellfun(@minus, coronal_cy, ...
    mat2cell(centroids_x, ones(1,T), ones(C,1)), ...
    'UniformOutput', 0);

%%

coronal_angles = cellfun(@atan2d, ccy_centered, ccx_centered, ...
    'UniformOutput',0);
% Correct for orientation of angles -- sort clock/anti-clockwise
orientedIndex = cell(T,C);
for t = 1:T
    for i = 1:C
        
        [~,orientedIndex{t,i}] = sort(coronal_angles{t,i});
        
    end
end

%%

diffX = cellfun(@(x,y) diff(x(y)), ccx_centered, orientedIndex, ...
    'UniformOutput',0);
diffY = cellfun(@(x,y) diff(x(y)), ccy_centered, orientedIndex, ...
    'UniformOutput',0);

coronal_perim = cellfun(@(x,y) nansum(sqrt(x.^2 + y.^2)), diffX, diffY);

perim = fits_twist.get_corrected_measurement(cells,coronal_perim,input);

%%

% embryoID = 1:5;
embryoID = 6:10;

for c = 1:3
    
    fOI = fits.get_embryoID(embryoID).get_cluster(c);
    fOI = fOI.find_near_fits(cellsOI,neighbor_defition);
    cellsOI = cells.get_embryoID(embryoID);
    
    for bin = 1:10
        
        fitsOI = fOI([fOI.bin] == bin);
        
        nearIDs = cat(1,fitsOI.nearIDs);
        num_near = cellfun(@(x) numel(x(~isnan(x))), nearIDs);
        
        % fitsOI.get_cluster(1)
        perim = fitsOI.get_corrected_measurement(cellsOI,coronal_perim,input);
        
        subplot(3,1,c)
        foo = cell(1,16);
        for i = 0:15
            foo{i+1} = nanmean(perim(num_near(:,3) == i,:),2);
        end
        
        [foo{cellfun(@isempty,foo)}] = deal(NaN);
        % throw out only 1-element parts
        [foo{cellfun(@(x) numel(x) < 2,foo)}] = deal(NaN);
        
        hold on
        subplot(3,1,c);
        medians(:,bin) = cellfun(@nanmedian,foo);
        plot(medians(:,bin),'Color',C(bin,:));
        xlim([0 10])
        
        %     errorbar(nanmedian(perim),nanstd(perim));
        
    end
end

