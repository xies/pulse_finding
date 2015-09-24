function [freq,center] = get_frequency(pulse)
% GET_FREQUENCY Estimate the frequency of pulsing in a given
% set of cells,
%
% Usage: [freq,center] = pulse.get_frequency;
%
% OUTPUT: frequency - a 1xNcell cell array of waiting time
%                     between pulses for each cell
%         center - the pairwise mean timing of pairs of
%                  consecutive pulses

% Return Fitted obj in cell arrays for each CellObj
fits = [pulse.fits]; cells = [pulse.cells];
cells = mat2cell(cells,1,ones(1,numel(cells)));

fits_incell = cellfun(@pulse.find_fits_from_cell, cells, ...
    'UniformOutput',0);

fits_center_incell = cell(1,numel(fits_incell));
for i = 1:numel(fits_incell)
    % Return each fit's center
    fits_incell{i} = fits_incell{i}.sort('center');
    fits_center_incell{i} = [fits_incell{i}.center];
end

% Use @diff to get waiting time intervals
freq = cellfun(@diff, fits_center_incell,'UniformOutput',0);
% "center" of a pulse pair (pair-wise means)
center = cellfun(@sort_pair_mean, fits_center_incell,'UniformOutput',0);

end