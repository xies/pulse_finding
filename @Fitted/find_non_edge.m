function f = find_non_edge(fits,cells)
%FIND_NON_EDGE Finds the fits that are not on the edge of the
% segmented embryo

embryos = unique([fits.embryoID]);

f = [];
for embryoID = embryos
    
    fitsOI = fits.get_embryoID(embryoID);
    cellsOI = cells.get_embryoID(embryoID).get_curated;
    cIDs = [cellsOI.cellID];
    
    tref = find( cellsOI(1).dev_time == 0);
    
    neighborhood = cat(2,cellsOI.identity_of_neighbors_all);
    neighborhood = neighborhood(tref,:);
    
    N = cellfun(@(x) numel(x(ismember(x,cIDs))), neighborhood);
    sIDs = [cellsOI(N > 3).stackID];
    
    f = [f fitsOI.get_stackID(sIDs)];
    
end

end % find_non_edge
