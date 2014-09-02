%% myosin_persistence.m

myos_diff = nanmin(corrected_myosin(:,5:end),[],2) - nanmin(corrected_myosin(:,1:4),[],2);
myos_diff_norm = myos_diff./nanmean(corrected_myosin(:,:),2);

fitsOI = fits.get_embryoID(1:5);

% myo_diff_norm = myos_diff_norm( ...
%     ismember([fits.embryoID],unique([fitsOI.embryoID])) );

behaviors = {'Ratcheted','Unratcheted','Unconstricting','N/A'};
l = [fitsOI.cluster_label];

% myo_diff_norm = myo_diff_norm(l < 4);
% l = l(l < 4);

labels = behaviors(l);

boxplot(myos_diff_norm,labels)
ylabel('Myosin persistence')
