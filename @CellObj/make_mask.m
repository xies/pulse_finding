function mask = make_mask(obj_array, frames, input)
%@CellObj.MAKE_MASK Make a BW mask of CellObjs using poly2mask.
%
% SYNOPSIS: mask = cells.make_mask(frames,input)
%
% INPUT:	cellobj - a array of cellobjs
%			frames - frames of interest
% 			input - input information (see LOAD_EDGE_SCRIPT)
% OUTPUT: 	mask - binary mask of size [X,Y,numel(frames)]
%
% See also: draw_measurement_on_cells, draw_measurement_on_cells_patch
%
% xies@mit August 2013

num_cells = numel(obj_array);
num_frames = numel(frames);

% if input.um_per_px not supplied, use 1
if ~isfield('input','um_per_px')
    um_per_px = 1;
else
    um_per_px = input.um_per_px;
end

X = input.X; Y = input.Y;
mask = zeros(Y,X,num_frames);

for i = 1:num_cells
    vx = obj_array(i).vertex_x; vy = obj_array(i).vertex_y;
    % etract frames of interest
    vx = vx(frames,:); vy = vy(frames,:);
    for t = 1:num_frames
        x = vx{t}/um_per_px; y = vy{t}/um_per_px;
        if all(~isnan(x))
            mask(:,:,t) = mask(:,:,t) + ...
                poly2mask(x,y,Y,X);
        end
    end
end
mask = logical(mask);

end % make_mask