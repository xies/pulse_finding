function adjust_dev_time(fits, old_tref, new_tref, dt)
%ADJUST_DEV_TIME Recalculate dev_time according to new reference point.
%
% USAGE: fits.adjust_dev_time( old_tref, old_tref, dt)
%
% Properties that will be adjusted:
%   fits.center
%   fits.dev_time

% Only one embryo
assert( numel(unique([fits.embryoID])) == 1, 'Fitted from only one embryo please.');

% Adjust temporal properties
for i = 1:numel(fits)
    
    this_fit = fits(i); % Passed by reference! No need to put back into array.
    
    this_fit.center = this_fit.center - (new_tref - old_tref)*dt;
    this_fit.dev_time = this_fit.dev_time ...
        - (new_tref - old_tref)*dt;
    
end

end