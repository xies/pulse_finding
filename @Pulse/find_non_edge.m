function f = find_non_edge(pulse)
%FIND_NON_EDGE Finds the fits that are not on the edge of the
% segmented embryo.
%
% USAGE: fits_non_edge = pulses.find_non_edge;

f = [];
for i = 1:numel(pulse);
    
    fitsOI = pulse(i).fits;
    cellsOI = pulse(i).cells.get_curated;
    
    for j = 1:numel(fitsOI)
        
        tref = fitsOI(j).center_frame;
        this_cellID = fitsOI(j).cellID;
        neighborhood = cellsOI([cellsOI.cellID] == this_cellID).identity_of_neighbors_all{tref};

        
        % Is # of curated neighbors
        N = ismember(neighborhood, [cellsOI.cellID]);
        N = numel(N(N > 0));
        
        % need all neighbors to be curated AND min of 3 neighbors
        if N > 3 && N == numel(neighborhood)
            f = [f fitsOI(j)];
        end
        
    end
    
end
    
end % find_non_edge
