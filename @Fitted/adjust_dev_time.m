function fits = adjust_dev_time(fits, old_tref, new_tref, dt)
%ADJUST_DEV_TIME Recalculate dev_time according to new reference point
% Properties that will be adjusted:
%   center
%   dev_time

% Only one embryo
if numel(unique([fits.embryoID])) > 1
    error('Fitted from only one embryo please.');
end

for i = 1:numel(fits)
    this_fit = fits(i);
    
    this_fit.center = this_fit.center - (new_tref - old_tref)*dt;
    this_fit.dev_time = this_fit.dev_time ...
        - (new_tref - old_tref)*dt;
    
    fits(i) = this_fit;
end
end