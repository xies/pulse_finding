function nearby_cells = get_nearby(obj_array,stackID,radius,reference_frame)
%@CellObj.GET_NEARBY Returns cells within radius r of a given
% cell.
%
% SYNOPSIS: nearby_cells = cells.get_nearby(stackID,radius,ref_frame)
%
% INPUT: stackID - the central cell
%        radius - the radius cutoff
%        ref_frame - the reference frame in the movie to use
%
% OUTPUT: nearby_cells
% xies@mit Feb 2014

central_cell = obj_array.get_stackID(stackID);
% Find all cells in the same embryo that's not the central cell
same_embryo = obj_array.get_embryoID( central_cell.embryoID);
same_embryo = same_embryo([same_embryo.stackID] ~= central_cell.stackID);

cx = [same_embryo.centroid_x]; cx = cx(reference_frame,:);
cy = [same_embryo.centroid_y]; cy = cy(reference_frame,:);

d = sqrt( ...
    (cx - central_cell.centroid_x(reference_frame)).^2 ...
    + (cy - central_cell.centroid_y(reference_frame)).^2);

nearby_cells = same_embryo( d <= radius );

end % get_nearby