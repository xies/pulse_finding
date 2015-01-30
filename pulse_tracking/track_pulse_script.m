%%TRACK_PULSE_SCRIPT Pipeline for

%%

mdf_file{6} = '~/Dropbox (MIT)/Adam/Tracked pulses/11-10-2012-3/11-10-2012-3_mimi.mdf'; embryoID(6) = 6;
% mdf_file{7} = '~/Dropbox (MIT)/Adam/Tracked pulses/cta pulses/twist_injection_022_mimi.mdf'; embryoID(7) = 7;
% mdf_file{8} = '~/Dropbox (MIT)/Adam/Tracked pulses/cta pulses/merge_t1-70.mdf'; embryoID(8) = 8;
mdf_file{9} = '~/Dropbox (MIT)/Adam/Tracked pulses/01-31-2012-1/merge_z6.mdf'; embryoID(9) = 9;

match_thresh = 1;

for i = 9
    
    % Load MDF into matrix
    mdf_mat = read_mdf(mdf_file{i});
    [tracks,cells_raw] = load_mdf_track(mdf_mat, embryo_stack, embryoID(i), ...
        input(i).t0, cells_raw);
    
    % Perform matching to fitted pulses
    
    pulse(i) = Pulse(tracks,mdf_file{i},fits_raw,fit_opts,cells_raw,input);
    pulse(i) = pulse(i).match_pulse(match_thresh);
    pulse(i) = pulse(i).categorize_mapping;
    
    pulse(i).embryoID = embryoID(i);
    display(['EmbryoID: ' num2str(i)]);
    display(mdf_file{i});
    display(pulse(i));
    
end

fits_curated = [pulse.fits];
cells_curated = [pulse.cells];
