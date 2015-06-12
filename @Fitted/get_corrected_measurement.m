function M = get_corrected_measurement(fits,c,meas,input)
fits = fits.align_fits(c,'measurement',meas);
fits = fits.resample_traces('measurement',[input.dt]);
M = cat(1,fits.corrected_measurement);
end