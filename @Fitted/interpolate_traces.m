function fits = interpolate_traces(fits,name,dt)
%interpolate_traces Uses INTERP1 to resample short traces
%
% [aligned_traces,aligned_time] = interpolate_traces(traces,embryoID,dt);
%
% xies@mit.edu Oct 2012

% concatenate traces
traces = cat(1,fits.(name));

% get dimensions
num_traces = size(traces,1);
embryoIDs = [fits.embryoID];

% Check to see if it's an cell array or not
if iscell(traces)
    display('Cannot do resampling on cell array measurements; returning original');
end

if numel(embryoIDs) ~= num_traces
    error('The number of traces and the number of embryoID must be the same.');
end

% find the aligned DT
aligned_dt = round(mean(dt)*100)/100;
% w = floor(T/2);

% Resample using the SIGNAL_PROCESSING TOOLBOX
for i = 1:num_traces
    
    if iscell(traces)
        fits(i).(['corrected_' name]) = traces(i,:);
        continue
    end
    
    l = fits(i).opt.left_margin;
    r = fits(i).opt.right_margin;
    
    % aligned_traces = zeros([num_traces, l + r - 3]);
    aligned_t = (- l : r )*aligned_dt;
    trace = traces(i,:);
    
    trace = trace( ~isnan(trace) );
    
    x = (-l:r)*dt( embryoIDs(i) );
    x = x( ~isnan(trace) );
    
    if numel(x) > 2
        % recenter
        x = interp1( x, trace, (-(l-2):r-2)*aligned_dt );
        if ~strcmpi(name,'myosin') && ~strcmpi(name,'measurement')
            fits(i).(['corrected_' name]) = x - nanmean(x);
        else
            fits(i).(['corrected_' name]) = x;
        end
    else
        fits(i).(['corrected_' name]) = ...
            nan(1, l+r-3 );
    end
    
    fits(i).corrected_time = aligned_t(3:end-2);
    
end

end % resample_traces