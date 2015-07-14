function cells_output = find_cells_with_fit(pulse,fits)
%FIND_CELL_WITH_FIT Return the CellObjs that contains the Fitted objects.
%
% USAGE: cells = pulses.find_cell_with_fit(fits);

embryoIDs = unique([fits.embryoID]);
assert( all(ismember(embryoIDs,[pulse.embryoID])),'Not all embryoIDs found in Pulse input');

cells_output = cell(1,numel(embryoIDs));

for i = 1:numel(embryoIDs)

    fitsOI = fits([fits.embryoID] == i);
    cells = [pulse(i).cells];
    cells = cells( ismember([cells.cellID],[fitsOI.cellID]) );

    cells_output{i} = cells;
    
end

cells_output = [cells_output{:}];

end