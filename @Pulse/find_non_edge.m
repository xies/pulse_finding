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
        neighborhood = cat(2,cellsOI.identity_of_neighbors_all);
        neighborhood = neighborhood(tref,:);
        
        N = cellfun(@(x) numel(x(ismember(x,[cellsOI.cellID]))), ...
            neighborhood);
        
%         foo = pulse(i).find_fits_from_cell(cellsOI(N > 3));
        if N([cellsOI.cellID] == fitsOI(j).cellID) > 3
            f = [f fitsOI(j)];
        end
        
    end
    
end
    
end % find_non_edge
