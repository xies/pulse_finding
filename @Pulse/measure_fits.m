function measure_fits(pulse)
%MEASURE_FITS Populates subsequence from CellObj data (EDGE data) into
% FITTED objects.
%
% USAGE: pulse.measure_fits;

fits = [pulse.fits];

pulse.align_fits('myosin');
pulse.align_fits('area');

aligned_area = cat(1,fits.area);
aligned_myosin = cat(1,fits.myosin);
% aligned_area_rate = cat(1,fits.area_rate);
% aligned_myosin_rate = cat(1,fits.myosin_rate);
% aligned_measurement = cat(1,fits.measurement);

% Mean-center pulses responses
aligned_area_norm = bsxfun(@minus,aligned_area,nanmean(aligned_area,2));
% aligned_myosin = bsxfun(@minus,aligned_myosin,nanmean(aligned_myosin,2));
% aligned_myosin = bsxfun(@rdivide,aligned_myosin,nanstd(aligned_myosin,[],2));
% aligned_measurement = bsxfun(@minus,aligned_measurement,nanmean(aligned_measurement,2));
pulse.assign_datafield(aligned_area_norm,'area_norm');
pulse.assign_datafield(aligned_myosin,'myosin');
% fits = assign_datafield(fits,aligned_measurement,'measurement');

% correlate for framerate differences
pulse.interpolate_traces('area_norm');
pulse.interpolate_traces('area');
pulse.interpolate_traces('myosin');

end