function W = get_pulse_transition_matrix(pulse,varargin)
%GET_PULSE_TRANSITION_MATRIX Constructs a matrix of the transition
% rates of different types of pulses (pooled across time).
%
% USAGE: W = pulse.get_pulse_transition_matrix

% fits = [pulse.fits];
if nargin < 2
    cells = [pulse.cells];
else
    cells = varargin{1};
end
% Filter fitted pulses by cells
fits = pulse.find_fits_from_cell( cells );

num_behavior = numel(unique([fits.cluster_label]));
% num_behavior = numel(unique([fits.bin]));
if num_behavior == 0
    error('No behavior labels are defined.')
end

% preallocate
W = zeros(num_behavior);

for i = 1:numel(fits)
    
    this_fit = fits(i);
    this_cell = pulse.find_cells_with_fit( this_fit );
    
    % fitID within this cell
    all_fits = pulse.find_fits_from_cell( this_cell ).sort('center');
    idx = find( [all_fits.fitID] == this_fit.fitID );
    
    % if this fit is not the last one
    if idx < numel(all_fits)
        next_label = all_fits( idx + 1 ).cluster_label;
%         next_label = all_fits( idx + 1 ).bin;
        
        W( this_fit.cluster_label, next_label ) = ...
            W( this_fit.cluster_label, next_label ) + 1;
%         W( this_fit.bin, next_label ) = ...
%             W( this_fit.bin, next_label ) + 1;
    end
end

end