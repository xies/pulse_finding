function [cx,cy,ct] = export_xyt( pulse, filename, traceback)
% EXPORT_XYT Exports the spatial and temporal coordinates of all pulses
% 
% NB: trackback - turn 'on' to use the earliest tracked
% cell centroid instead of instantaneous cell locations
%
% USAGE: [cx,cy,ct] = pulse.export_xyt(filename, traceback)


if nargin < 4, traceback = 'off'; end

fits = [pulse.fits];
cells = [pulse.cells];

cx = zeros(1,numel(fits));
cy = zeros(1,numel(fits));

for i = 1:numel(fits)
    
    this_fit = fits(i);
    
    x = cells.get_fitID( this_fit.fitID ).centroid_x;
    y = cells.get_fitID( this_fit.fitID ).centroid_y;
    
    if strcmpi(traceback,'on')
        
        cx(i) = x(find_earlierst_nonan(x));
        cy(i) = y(find_earlierst_nonan(y));
        
    else
        
        cframe = this_fit.center_frame;
        
        cx(i) = x( cframe );
        cy(i) = y( cframe );
        
        if isnan(cx(i))
            I = find_nearest_nonan( x, cframe );
            cx(i) = x(I);
            cy(i) = y(I);
            if isnan(cx(i)), keyboard; end
        end
    end
    
end

ct = [fits.center];
l = [fits.cluster_label];
fIDs = [fits.fitID];

[ct,order] = sort(ct,'ascend');

if ~isempty(filename)
    M = cat(1,fIDs,cx(order),cy(order),ct,l)';
    csvwrite(filename,M);
end

end % export_xyt