function myosin_persistence = get_myosin_persistence(fits)
% GET_MYOSIN_PERSISTENCE

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