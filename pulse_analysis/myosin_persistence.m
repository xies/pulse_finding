%% myosin_persistence.m

fitsOI = fits.get_embryoID(6:10);

myo_persis = get_myosin_persistence(fitsOI);

behaviors = {'Ratcheted','Unratcheted','Unconstricting','N/A'};
l = [fitsOI.cluster_label];

% myo_diff_norm = myo_diff_norm(l < 4);
% l = l(l < 4);

labels = behaviors(l);

boxplot(myo_persis,labels);
% distributionPlot(myos_diff_norm,'groups',labels)
ylabel('Myosin persistence')

%% Correlate ratcheting w/ persistence or # near pulses

% Define neighborhood
neighbor_definition.temporal.def = @(time_diff,tau) abs(time_diff) < tau & time_diff > 0;
neighbor_definition.temporal.windows = time_windows;
neighbor_definition.spatial.def = 'identity';
