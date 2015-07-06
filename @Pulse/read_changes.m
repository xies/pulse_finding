function pulse = read_changes( pulse, changes )
% READ_CHANGES Make edits to track/fit given recorded changes
% For use from MATCH_VIEWER

if isfield(changes,'fitsMadeFromTrack')
    tracks = [changes.fitsMadeFromTrack.tracks];
    %                 trackID = [trackIDs.trackID];
    for i = 1:numel(tracks)
        opt = pulse.fit_opt;
        trackID = pulse.find_nearest_object('track', ...
            tracks(i).cx,tracks(i).cy,tracks(i).ct).trackID;
        pulse = pulse.createFitFromTrack(trackID,opt);
    end
end

if isfield(changes,'tracksMadeFromFit')
    fits = [changes.tracksMadeFromFit.fits];
    %                 fitIDs = [fitIDs.fitID];
    for i = 1:numel(fits)
        fitID = pulse.find_nearest_object('fit', ...
            fits(i).cx,fits(i).cy,fits(i).ct);
        fitID = fitID.fitID;
        pulse = pulse.createTrackFromFit(fitID);
    end
end
if isfield(changes,'fitsRemoved')
    fits = [changes.fitsRemoved.fits];
    for i = 1:numel(fits)
        fitID = pulse.find_nearest_object('fit', ...
            fits(i).cx,fits(i).cy,fits(i).ct);
        if isempty(fitID)
            warning(['Can''t find fit at x = ' ...
                num2str(tracks(i).cx) ...
                ', y = ' num2str(tracks(i).cy) ...
                ', t = ' num2str(tracks(i).ct) ...
                ]);
        else
            fitID = fitID.fitID;
            pulse = removePulse(pulse,'fit',fitID);
        end
    end
end
if isfield(changes,'tracksRemoved')
    tracks = [changes.tracksRemoved.tracks];
    for i = 1:numel(tracks)
        trackID = pulse.find_nearest_object('track', ...
            tracks(i).cx,tracks(i).cy,tracks(i).ct);
        if isempty(trackID)
            warning(['Can''t find track at x = ' ...
                num2str(tracks(i).cx) ...
                ', y = ' num2str(tracks(i).cy) ...
                ', t = ' num2str(tracks(i).ct) ...
                ]);
        else
            trackID = trackID.trackID;
            pulse = removePulse(pulse,'track',trackID);
        end
    end
end
if isfield(changes,'reassignedTrackFit')
    changes = [changes.reassignedTrackFit];
    for i = 1:numel(changes)
        track = changes(i).track;
        fit = changes(i).fit;
        trackID = pulse.find_nearest_object('track', ...
            track.cx,track.cy,track.ct).trackID;
        fitID = find_nearest_object(pulse,'fit', ...
            fit.cx,fit.cy,fit.ct).fitID;
        pulse = reassignFit(pulse,fitID,trackID);
    end
end

save( [fileparts(pulse.tracks_mdf_file), '/', 'pulse_curated.mat'], 'pulse');

end % read_changes