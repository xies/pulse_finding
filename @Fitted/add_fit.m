function [obj_array,errorflag] = add_fit(obj_array,new_fit)
%ADD_FIT tries to append a new_fit to a FITS array. Fails if
% any fit in the array is equal to the new_fit.
%
% See also: FITTED.eq
errorflag = 0;
if any(obj_array == new_fit)
    disp('Cannot create new fit: Fit already exists.');
    beep
    errorflag = 1;
    return
end

obj_array = [obj_array new_fit];

end % add_fit