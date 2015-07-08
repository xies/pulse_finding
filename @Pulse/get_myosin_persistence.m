function myosin_persistence = get_myosin_persistence(pulse)
% GET_MYOSIN_PERSISTENCE
% Defined as the normalized difference between the pre-peak minimum myosin
% intensity and post-peak minimum myosin intensity during a pulse.
% 
% USAGE: MP = pulse.get_myosin_persistence;

fits = [pulse.fits];

corrected_myosin = cat(1,fits.corrected_myosin);
if size(corrected_myosin,1) ~= numel(fits);
    error('Not all fits have myosins associated with them.');
end

l = fits(1).opt.left_margin + 1;

myosin_persistence = nanmin(corrected_myosin(:,l:end),[],2) ...
    - nanmin(corrected_myosin(:,1:l),[],2);
myosin_persistence = ...
    myosin_persistence./nanmean(corrected_myosin(:,:),2);

end