function pulse = match_pulse(pulse,threshold)
%MATCH_PULSE Match FITTED objects to TRACK objects by generating a
% MatchTrackFit object, which is a wrapper for interacting with two
% container.Map objects (hashmaps).
%
% USAGE: pulse = match_pulse(pulse, threshold)

if isempty(pulse.map) || ...
        (~isempty(pulse.match_thresh) && pulse.match_thresh ~= threshold)
    pulse.categories = [];
    pulse.changes = [];
end

nbm = MatchTrackFit(pulse.tracks,pulse.fits,threshold);
pulse.map = nbm;
pulse.match_thresh = threshold;

end % match_pulse
