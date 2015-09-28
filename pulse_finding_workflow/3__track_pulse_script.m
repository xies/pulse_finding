%%TRACK_PULSE_SCRIPT Pipeline for

%% Load tracks into Pulse

mdf_file{1} = '~/Dropbox (MIT)/cta/raw_data/05-22-2015-5/merged_z7.mdf'; embryoID(1) = 1;
mdf_file{2} = '~/Dropbox (MIT)/cta/raw_data/05-22-2015-4/merged_z6.mdf'; embryoID(2) = 2;

mdf_file{3} = '~/Dropbox (MIT)/cta/raw_data/11-10-2012-3/11-10-2012-3_mimi.mdf'; embryoID(3) = 3;
mdf_file{4} = '~/Dropbox (MIT)/cta/raw_data/01-31-2013-1/merge_z6.mdf'; embryoID(4) = 4;

mdf_file{5} = '~/Dropbox (MIT)/cta/raw_data/08-23-2015-1/merged_z8.mdf'; embryoID(5) = 5;

match_thresh = 1;

for i = 4
    
    eID = embryoID(i);
    cellsTmp = cells_raw([cells_raw.embryoID] == eID);
    
    % Load MDF into matrix
    mdf_mat = read_mdf(mdf_file{i});
    tracks = load_mdf_track(mdf_mat, embryo_stack(eID), ...
        input(i).t0, cellsTmp);
    
    % Perform matching to fitted pulses
    pulse(i) = pulse(i).match_tracks_to_fits(tracks,mdf_file{i},match_thresh);
%     pulse(i) = pulse(i).match_pulse(match_thresh);
    pulse(i) = pulse(i).categorize_mapping;
    
    display(['EmbryoID: ' num2str(eID)]);
    display(pulse(i));
    
end

%% Curate pulses

match_viewer(pulse(4),embryo_stack(4) );
