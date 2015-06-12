
function [adj,nodes] = get_pulse_transition_graph(cells,fits)
%GET_PULSE_TRANSITION_GRAPH Constructs a matrix of the transition
% rates of different types of pulses (pooled across time).

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
    % fitID within this cell
    this_cellfitIDs = cells.get_stackID( this_fit.stackID).fitID;
    idx = find( this_cellfitIDs == this_fit.fitID );
    % if this fit is not the last one
    if idx < numel(this_cellfitIDs)
        next_fitID = this_cellfitIDs( idx + 1 );
        next_label = fits.get_fitID( next_fitID ).cluster_label;
        
        adj( this_fit.cluster_label, next_label + num_behavior ) = ...
            adj( this_fit.cluster_label, next_label + num_behavior ) + 1;
    end
end

end