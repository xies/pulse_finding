function fits_new = sort(fits,field)
% SORT pulses by a given field, Default = amplitude (ascending)
% Will not do in-place sort.
%
% USAGE: fits_sorted = sort(fits,field2sort);
fits_new = fits.copy;
if nargin < 2, field = 'amplitude'; end
% if strcmpi(field,'cluster_weight')
%     U = cat
%     tobeSorted = 
% else
    tobeSorted = cat(1,fits.(field));
% end
[~,order] = sort( nanmax( tobeSorted,[], 2),'ascend');
fits_new = fits_new(order);
end % sort