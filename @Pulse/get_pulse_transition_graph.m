function [adj,nodes] = get_pulse_transition_graph(pulse,varargin)
%GET_PULSE_TRANSITION_GRAPH Constructs a matrix of the transition
% rates of different types of pulses (pooled across time).

if nargin < 2
    cells = [pulse.cells];
else
    cells = varargin{1};
end

fits = pulse.find_fits_from_cell( cells );
num_behavior = numel(unique([fits.cluster_label])); 

% construct the nodes
nodes = cat(1, ...
    cat(2,ones(num_behavior,1),(1:num_behavior)'), ...
    cat(2,2*ones(num_behavior,1),(1:num_behavior)') ...
    );

% preallocate the adjacency matrix
adj = zeros(2*num_behavior);
for i = 1:numel(fits)
    
    this_fit = fits(i);
    this_cell = pulse.find_cells_with_fit( this_fit );
    
    % fitID within this cell
    all_fits = pulse.find_fits_from_cell( this_cell ).sort('center');
    idx = find( [all_fits.fitID] == this_fit.fitID );
    
    % if this fit is not the last one
    if idx < numel(all_fits)
        next_label = all_fits(idx + 1).cluster_label;
        % +1 for transition matrix
        adj( this_fit.cluster_label, next_label + num_behavior ) = ...
            adj( this_fit.cluster_label, next_label + num_behavior ) + 1;
    end
end

end