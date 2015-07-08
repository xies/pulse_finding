%%TRACK_PULSE_SCRIPT Pipeline for

%%

mdf_file{1} = '~/Dropbox (MIT)/cta/05-22-2015-5/merged_z7_2-5std.mdf'; embryoID(1) = 1;
mdf_file{2} = '~/Dropbox (MIT)/cta/05-22-2015-4/merged_z6.mdf'; embryoID(2) = 2;

match_thresh = 1;

for i = 1
    
    % Load MDF into matrix
    mdf_mat = read_mdf(mdf_file{i});
    [tracks,cells_raw] = load_mdf_track(mdf_mat, embryo_stack, embryoID(i), ...
        input(i).t0, cells_raw);
    
    % Perform matching to fitted pulses
    
    pulse(i) = Pulse(tracks,mdf_file{i},fits_raw,fit_opts,cells_raw,input(i));
    pulse(i) = pulse(i).match_pulse(match_thresh);
    pulse(i) = pulse(i).categorize_mapping;
    
    pulse(i).embryoID = embryoID(i);
    display(['EmbryoID: ' num2str(i)]);
    display(mdf_file{i});
    display(pulse(i));
    
end

fits_curated = [pulse.fits];
cells_curated = [pulse.cells];