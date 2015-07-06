
function pulse = removePulse(pulse,type,pulseID)
%@Pulse.removePulse Remove pulse from track-fit mapping, as well as
% the respective pulse object array.
%
% USAGE: pulse = removePulse(pulse,'fit',fitID);
%
% xies@mit.edu
new_nbm = pulse.map.removeElement(pulseID,type);
pulse.map = new_nbm;
c = pulse.cells;

% Remove pulse from stack
switch type
    case 'fit'
        indices = ismember([pulse.fitsOI_ID], pulseID);
        if ~any(indices)
            display('Cannot remove FITTED: given fitID does not exist.');
            return
        end
        % get a copy of to be deleted FITTED and record its
        % centroid + developmental time
        
        fit = pulse.fits.get_fitID( pulseID );
        c = c.get_stackID( fit.stackID );
        
        [cx,cy,ct] = fit.get_xyt(c);
        
        % Record b/f removal -- instead of its ID, use [cx,cy,ct]
        if isfield(pulse.changes,'fitsRemoved')
            pulse.changes.fitsRemoved.fits(end+1).cx = cx;
            pulse.changes.fitsRemoved.fits(end).cy = cy;
            pulse.changes.fitsRemoved.fits(end).ct = ct;
        else
            pulse.changes.fitsRemoved.fits.cx = cx;
            pulse.changes.fitsRemoved.fits.cy = cy;
            pulse.changes.fitsRemoved.fits.ct = ct;
        end
        
        % Remove from cell obj
        stackID = [pulse.fits.get_fitID(pulseID).stackID];
        pulse.cells( [pulse.cells.stackID] == stackID ) = ...
            pulse.cells( [pulse.cells.stackID] == stackID ).removeFit(pulseID);
        
        % Remove from fits stack
        pulse.fits = pulse.fits.removeFit( pulseID );
        pulse.fitsOI_ID(pulse.fitsOI_ID == pulseID) = [];
        
        display(['Deleting fitID: ' num2str(pulseID)]);
        
    case 'track'
        indices = ismember([pulse.tracks.trackID], pulseID);
        if ~any(indices)
            display('Cannot remove TRACK: given trackID does not exist.');
            return
        end
        
        % get a copy of to be deleted TRACK and record its
        % centroid + developmental time
        track = pulse.tracks.get_trackID( pulseID );
        c = c.get_stackID( track.stackID );
        
        [cx,cy,ct] = track.get_xyt(c);
        
        % Record b/f removal - [cx,cy,ct]
        if isfield(pulse.changes,'tracksRemoved')
            pulse.changes.tracksRemoved.tracks(end+1).cx = cx;
            pulse.changes.tracksRemoved.tracks(end).cy = cy;
            pulse.changes.tracksRemoved.tracks(end).ct = ct;
        else
            pulse.changes.tracksRemoved.tracks.cx = cx;
            pulse.changes.tracksRemoved.tracks.cy = cy;
            pulse.changes.tracksRemoved.tracks.ct = ct;
        end
        
        % Remove from CellObj
        stackID = [pulse.tracks(indices).stackID];
        for i = 1:numel(stackID)
            pulse.cells( [pulse.cells.stackID] == stackID(i)) = ...
                pulse.cells( [pulse.cells.stackID] == stackID(i)).removeTrack(pulseID);
        end
        % Remove from track stack
        pulse.tracks( indices ) = [];
        display(['Deleting trackID: ' num2str(pulseID)]);
        
    otherwise
        error('Invalid type: expecting TRACK or FIT.')
end

% Redo categorizing
pulse = pulse.match_pulse(pulse.match_thresh);
pulse = pulse.categorize_mapping;

end %removePulse