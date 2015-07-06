function corona_measurement = get_corona_measurement(cells,measurement)
%GET_CORONAL_MEASUREMENT - return the given measurement
% over the +1 neighbors (corona) on a cell array

%             [num_frames,num_cells] = size(measurement);
num_frames = numel(cells(1).dev_time);
embryoID = unique([cells.embryoID]);

corona_measurement = cell(num_frames,numel(cells));
cellID_padding = 0;

for eID = embryoID
    
    cells_in_embryo = cells.get_embryoID(eID);
    num_cells = numel(cells_in_embryo);
    this_m = cell(num_frames,num_cells);
    
    for t = 1:num_frames
        
        for i = 1:num_cells
            
            IDs = cells_in_embryo(i).identity_of_neighbors_all{t};
            % pad cellID so we get right indexing into the
            % concatenated measurements matrix
            IDs = IDs + cellID_padding;
            
            if ~isnan(IDs)
                IDs = IDs(IDs > 0);
                %                             this_m(t,i) = nansum(measurement(t,IDs));
                this_m{t,i} = measurement(t,IDs);
            end
            
        end
    end
    
    corona_measurement(:,[cells.embryoID] == eID) ...
        = this_m;
    cellID_padding = cellID_padding + num_cells;
    
end
end % get_coronal_measurement