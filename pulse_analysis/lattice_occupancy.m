function R = lattice_occupancy(cells)
%LATTICE_OCCUPANCY


N = numel(cells);
T = size( cells(1).area,1);

R = cell(N,T);

for i = 1:N  % iterate through all cells
    for t = 1:T % iterate through each time slice
        
        this_cell = cells(i); 
%         dframe = fix( timewindow/nanmean(diff(this_cell.dev_time)) );
        
        if ~isnan(this_cell.dev_time(t))
            
            same_embryo = cells.get_embryoID( this_cell.embryoID );  % find all the cells in the same embryo
            % not the same cell (self)
            same_embryo = same_embryo( [same_embryo.cellID] ~= this_cell.cellID); % exclude current "this_cell"
            
            cx = this_cell.centroid_x(t);
            cy = this_cell.centroid_y(t);
            
            nx = cat(2,same_embryo.centroid_x);  % put position of all other cells in matrix, dimension (row- time x col-cells).
            ny = cat(2,same_embryo.centroid_y);
            nx = nx( t,: );                 % take all cells for certain time point by grabbing row
            ny = ny( t,: );
            
            lx = nx - cx;                   % subtract x-position of this_cell
            ly = ny - cy;
            
            foo = sqrt( lx.^2 + ly.^2 );    % distance between cell of interest and all other cells 
            R{i,t} = cat(1, 0, foo(:));     % add a zero before the distance vector 'foo'
            
            % we have a fector of distances to neighbors, but how do we
            % keep track of what cell that corresponds to?
            
        end
        
    end
    display(['Done with cell #' num2str(i)]);
end