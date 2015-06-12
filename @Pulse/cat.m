function pulse = cat(pulse1,pulse2)
%CAT - overloaded concatenation method for putting together two
% Pulse objects. Will look through embryoIDs to ensure there
% are no colliding ID numbers.

pulse = pulse1;
embryoIDs1 = [pulse1.embryoID]; embryoIDs2 = [pulse2.embryoID];

% find colliding/not colliding embryoIDs
conflicting_embryoIDs = union( embryoIDs1, embryoIDs2 );

if isempty( conflicting_embryoIDs)
    pulse = [pulse pulse2];
else
    
    for i = 1:numel(conflicting_embryoIDs)
        
        this_embryo_pulse = [pulse2.embryoIDs];
        % new embryoID is +1 to max embryoIDs
        new_embryoID = max([embryoIDs1,embryoIDs2]) + 1;
        % reindex fitIDs
        this_embryo_pulse.fits = ...
            this_embryo_pulse.fits.reindex_fitID(new_embryoID);
        % reindex trackIDs
        this_embryo_pulse.tracks = ...
            this_embryo_pulse.tracks.reindex_trackID(new_embryoID);
        
        pulse = [pulse this_embryo_pulse];
        
    end
    
end

end % cat