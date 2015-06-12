function fits = sort(fits,field)
% SORT pulses by a given field, Default = amplitude (ascending)
%
% USAGE: fits_sorted = sort(fits,field2sort);
if nargin < 2, field = 'amplitude'; end
[~,order] = sort( nanmax( cat(1,fits.(field)),[], 2));
fits = fits(order);
end % sort
