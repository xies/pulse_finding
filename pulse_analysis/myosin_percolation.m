fitsOI = fits_wt;

tot_num_cell = num_cells([fitsOI.embryoID]);
M = fitsOI.get_corrected_measurement(cells,myosin_perc,input);
M = bsxfun(@rdivide,M,tot_num_cell');

%%

% color = hsv(3);
for i = 1:3
    
    subplot(1,3,i);
    x = fits(1).corrected_time;
    y = M([fitsOI.cluster_label] == i,:);
    plot(x,nanmean(y))
    
    ylim([0 1])
    
end

%%

color = hsv(3);
for i = 1:3
    
    subplot(3,1,i);
    y = M([fitsOI.cluster_label] == i,:);
    y = mean(y,2);
    
    I = [fitsOI.cluster_label] == i;
    
    scatter([fitsOI(I).center],y);
    ylim([0 1]), xlim([-300 300])
    
end
