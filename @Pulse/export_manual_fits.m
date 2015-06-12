function export_manual_fits(pulse)
%EXPORT_MANUAL_PULSES
% Writes down the manual fit parameters for fits created from
% tracks. Saveas as manual_fit.csv
changes = pulse.changes;
if isfield(changes,'fitsMadeFromTrack')
    num_changes = numel(changes.fitsMadeFromTrack);
else
    num_changes = 0;
end

mat2write = nan(num_changes,6);

for i = 1:num_changes
    this_change = changes.fitsMadeFromTrack(i);
    
    %                 trackID = this_change.trackID;
    %                 [cx,cy,ct] = STlocation(pulse,pulse.tracks.get_trackID(this_change.trackID) );
    fitID = pulse.find_nearest_object('fit',...
        this_change.fits.cx,this_change.fits.cy,this_change.fits.ct).fitID;
    
    this_fit = pulse.fits.get_fitID(fitID);
    if ~isempty(this_fit)
        params = [this_fit.amplitude this_fit.center this_fit.width];
        mat2write(i,1) = this_change.tracks.cx;
        mat2write(i,2) = this_change.tracks.cy;
        mat2write(i,3) = this_change.tracks.ct;
        mat2write(i,4:6) = params;
    end
end

if num_changes > 0
    csvwrite( [fileparts(pulse.tracks_mdf_file), '/', 'manual_fits.csv'], ...
        mat2write );
end

end % export_manual_fits
