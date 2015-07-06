function [fits] = assign_datafield(fits,data,name)
%assign_datafield Assign a matrix (Nc x p) into the fields of an array of
% Fitted. The size of matrix and number of fits must match.
%
% USAGE: fits.assign_datafield(data,name)
%
if size(data,1) ~= numel(fits)
    error('Data size must be the same as the number of FITTED objects.');
end
for i = 1:numel(fits)
    fits(i).(name) = data(i,:);
end
end % assign_datafield
