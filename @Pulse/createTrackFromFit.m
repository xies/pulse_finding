
function pulse = createTrackFromFit(pulse,fitID)
%@Pulse.createTrackFromFit Convert a fitted pulse into an
%'artificial track'. Useful when dealing with the 'add'
% category.
%
% USAGE: pulse = pulse.createTrackFromFit(fitID);

fit = pulse.fits.get_fitID(fitID);
if isempty(fit), display('No FIT found with fitID.'); return; end
if isfield(pulse.changes,'tracksMadeFromFit')
    changes = [pulse.changes.tracksMadeFromFit.fits];
    already_used = zeros(1,numel(changes));
    for i = 1:numel(changes)
        already_used(i) = pulse.find_nearest_object('fit', ...
            changes(i).cx,changes(i).cy,changes(i).ct).fitID;
    end
    if any( fitID == already_used)
        display(['Fit #' num2str(fitID) ' already used to add a track.']);
        return
    end
end

display(['Creating track from fitID ' num2str(fitID)]);

% Add to tracks stack
this_track.embryoID = fit.embryoID;
this_track.cellID = fit.cellID;
this_track.stackID = fit.stackID;
this_track.dev_frame = fit.width_frames;
this_track.embryoID = fit.embryoID;
this_track.dev_time = ensure_row(fit.dev_time);

tracks = add_track(pulse.tracks,this_track);

% Rematch the track/fit mappings
pulse_new = pulse;
pulse_new.tracks = tracks;
pulse_new = pulse_new.match_pulse(pulse.match_thresh); % Redo match
pulse_new = pulse_new.categorize_mapping;

pulse = pulse_new;

pulse.cells([pulse.cells.stackID] == fit.stackID) = ...
    pulse.cells([pulse.cells.stackID] == fit.stackID ).addTrack( pulse.tracks(end).trackID);

% Record changes
c = pulse.cells.get_stackID(fit.stackID);
[cx,cy,ct] = fit.get_xyt(c);
if isnan(cx), keyboard; end

this_change.tracks.cx = cx;
this_change.tracks.cy = cy;
this_change.tracks.ct = ct;
this_change.fits.cx = cx;
this_change.fits.cy = cy;
this_change.fits.ct = ct;

if isfield(pulse.changes,'tracksMadeFromFit')
    pulse.changes.tracksMadeFromFit(end+1) = this_change;
else
    pulse.changes.tracksMadeFromFit = this_change;
end

end % createTrackFromFit