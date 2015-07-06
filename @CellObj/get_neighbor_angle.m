function angles = get_neighbor_angle(cellx,celly,frame)
%GET_NEIGHBOR_ANGLE
% Given a pair of cells, get the angle between their centroids (y wrt x,
% along x-axis).
%
% Will flip all angles to lie within quadrant I. (i.e. use the
% leftmost cell as a reference). Will return degrees.
%
% Optionally, will return a angle at a specific frame, if
% supplied
%
% SYNOPSIS: angles = get_neighbor_angle(cell1,cell2);
%           angle = get_neighbor_angle(cell1,cell2,frame);
%
% xies@mit March 2012

% check that they're in the same embryo
if cellx.embryoID ~= celly.embryoID, error('Cells not in same embryo!'); end

angles = rad2deg(rad_flip_quadrant(atan2(...
    bsxfun(@minus,celly.centroid_y,cellx.centroid_y), ...
    bsxfun(@minus,celly.centroid_x,cellx.centroid_x) ) ) );

if nargin > 2, angles = angles(frame); end
if angles < -90 || angles > 90, keyboard; end

end