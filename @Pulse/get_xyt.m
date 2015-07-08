function [cx,cy,ct] = get_xyt(fit,cell)
validateattributes(fit,{'Fitted'},{'scalar'});
validateattributes(cell,{'CellObj'},{'scalar'});

cframe = findnearest(nanmean(fit.dev_time),cell.dev_time);
if numel(cframe) > 1, cframe = cframe(1); end

cx = repnan(cell.centroid_x);
cx = cx(cframe);

cy = repnan(cell.centroid_y);
cy = cy(cframe);

ct = mean(fit.dev_time);

end