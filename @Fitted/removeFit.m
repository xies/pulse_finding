function obj_array = removeFit(obj_array,fitID)
%removeFit Removes a fit of a given fitID from an array
% USAGE: obj_array = obj_array.removeFit(fitID)
obj_array([obj_array.fitID] == fitID) = [];
end % removeFit