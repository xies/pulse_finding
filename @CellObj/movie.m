function F = movie(cells,stackID,embryo_stack)
% MOVIE - make a movie of the cell
% USAGE: F = movie(cells,stackID,embryo_stack)
% xies@mit.edu

% Extract cell
this_cell = cells.get_stackID(stackID);
% Get its embryo
this_embryo = embryo_stack(this_cell.embryoID);
% Get IDs
h.cellID = this_cell.cellID;
h.input = this_embryo.input;

% Get vertices
h.vx = this_embryo.vertex_x;
h.vy = this_embryo.vertex_y;

h.border = 'on';
h.frames2load = find(~isnan(this_cell.dev_time));

% Check that there are the measurements you're looking for...
if ~isempty(this_cell.myosin_sm)
    h.channels = {'Membranes','Myosin','Membranes'};
else
    h.channels = {'Membranes','Rho Kinase thresholded'};
end

if this_cell.flag_fitted
    h.measurement = nan(numel(h.frames2load),3);
    I = find( ismember(this_cell.fit_time, this_cell.dev_time) );
    %                 foo = sum(this_cell.fit_gausses,2);
    foo = bsxfun(@rdivide, this_cell.fit_gausses, max(this_cell.fit_gausses) );
    foo = sum(foo,2);
    h.measurement(I,:) = foo(:,ones(1,3));
end

F = make_cell_img(h);

end % movie