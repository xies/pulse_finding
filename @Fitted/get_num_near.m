function num_near = get_num_near(fits,cells,neighbor_definition,window)
%GET_NUM_NEAR Returns the number of nearby fits given a CellObj array and
% the spatiotemporal neighbor defition.
%
% USAGE: num_near = fits.get_num_near(cells, neighbor_def, window);

fits = fits.find_near_fits(cells,neighbor_definition);

nearIDs = cat(1,fits.nearIDs);

% Convert to number of pulses
num_near = cellfun(@(x) numel(x(~isnan(x))), nearIDs);
num_near = num_near(:,window);

end
