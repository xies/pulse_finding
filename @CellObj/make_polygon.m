
function [x,y] = make_polygon(obj_array,t,input,filename)
%@CellObj.make_polygon Return a polygon object containing all
% cells within a cellobj array. Most useful for inter-connected
% clusters of cells.
%
% SYNOPSIS: [x,y,flag] = make_polygon(obj_array,t,filename)
%
% INPUT: cellobj_array - has to be from the same embryo
%        t - the frame of interest
%        input- either just embyro-specific input or an array
%               of input from EDGE_LOAD_SCRIPT
%        filename - (optional) For saving into a CSV file
% OUTPUT: x,y - coordinates of union of cells
%         flag - if there is hole, will raise flag
%
% See also: POLYBOOL
%
% xies@mit.edu Oct 2013

if numel( unique([obj_array.embryoID]) ) > 1
    error('Cells have to be from the same embryo!');
end

if numel(input) > 1
    input = input(obj_array(1).embryoID);
end
vx = cat(2,obj_array.vertex_x);
vy = cat(2,obj_array.vertex_y);

vx = cellfun(@(x) x*input.um_per_px,vx(t,:),'UniformOutput',0);
vy = cellfun(@(x) x*input.um_per_px,vy(t,:),'UniformOutput',0);

%             % y-direction is off
%             vy = cellfun(@(x) input.Y-x,vy,'UniformOutput',0);

warning('off','map:vectorsToGPC:noExternalContours');
x = vx{1}; y = vy{1};
for i = 2:numel(vx)
    [x,y] = polybool('union',x,y,vx{i},vy{i});
end
[x,y] = poly2cw(x,y);
x = x(1:end-1); y = y(1:end-1);

if nargin > 3,
    csvwrite(filename, cat(2,x,y));
end

end % make_polygon
