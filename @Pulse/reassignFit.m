function pulse = reassignFit(pulse,fitID,newTrackID)
%@Pulse.reassginFit Reassign FIT to TRACK: works only if neither
% FIT and TRACK have been assigned previously.
nbm = pulse.map.reassign(newTrackID,fitID);
pulse.map = nbm;
pulse = pulse.categorize_mapping;

this_change.trackID = newTrackID;
this_change.fitID = fitID;
if isfield(pulse.changes,'reassignedTrackFit')
    pulse.changes.reassignedTrackFit = ...
        [pulse.changes.reassignedTrackFit, this_change];
else
    pulse.changes.reassignedTrackFit = this_change;
end

end % reassignFit