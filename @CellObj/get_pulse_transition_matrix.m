function W = get_pulse_transition_matrix(cells,fits)
%GET_PULSE_TRANSITION_MATRIX Constructs a matrix of the transition
% rates of different types of pulses (pooled across time).
%
%

num_behavior = numel(unique([fits.cluster_label]));
% Filter fitted pulses by cells stack
fits = fits.get_stackID([cells.stackID]);

% preallocate
W = zeros(num_behavior);

for i = 1:numel(fits)
    
    this_fit = fits(i);
    % fitID within this cell
    this_cellfitIDs = cells.get_stackID( this_fit.stackID).fitID;
    idx = find( this_cellfitIDs == this_fit.fitID );
    % if this fit is not the last one
    if idx < numel(this_cellfitIDs)
        next_fitID = this_cellfitIDs( idx + 1 );
        next_label = fits.get_fitID( next_fitID ).cluster_label;
        
        W( this_fit.cluster_label, next_label ) = ...
            W( this_fit.cluster_label, next_label ) + 1;
    end
end

end