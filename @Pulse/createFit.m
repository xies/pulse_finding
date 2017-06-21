function pulse = createFit(pulse,cellID,opt)
% Manually create a Fitted from scratch.

c = pulse.cells(cellID);

params = manual_fit( ...
    [100 20],c);

% Construct a new FITTED object from parameters
new_fit = Fitted( c, params, pulse.next_fitID, opt);
pulse.next_fitID = pulse.next_fitID + 1; %increment fitIDs

% Add into stack
[fits,errorflag] = pulse.fits.add_fit(new_fit);
if errorflag, return; end
fits(end).manually_added = 1;

pulse.fitsOI_ID = [pulse.fitsOI_ID fits(end).fitID];
pulse.fits = fits;

% Redo match/categorizing
pulse = pulse.match_pulse(pulse.match_thresh);
pulse = pulse.categorize_mapping;

% Update cell tracklist
% pulse.cells( [pulse.cells.stackID] == track.stackID ) = ...
%     pulse.cells( [pulse.cells.stackID] == track.stackID).addFit( pulse.fits(end));
c.addFit(pulse.fits(end));

% % Record changes
% %             this_change.fitID = pulse.fits(end).fitID;
% %             this_change.trackID = trackID;
% [cx,cy,ct] = pulse.get_xyt(track);
% this_change.fits.cx = cx;
% this_change.fits.cy = cy;
% this_change.fits.ct = ct;
% this_change.tracks.cx = cx;
% this_change.tracks.cy = cy;
% this_change.tracks.ct = ct;
% 
% if isfield(pulse.changes,'fitsMadeFromTrack')
%     pulse.changes.fitsMadeFromTrack(end+1) = this_change;
% else
%     pulse.changes.fitsMadeFromTrack = this_change;
% end

end % createFit