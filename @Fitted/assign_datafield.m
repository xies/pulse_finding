function [fits] = assign_datafield(fits,data,name)
%assign_datafield Assign a matrix into the fields of fitted objects
% USAGE: fits = assign_datafield(fits,data,name)
if size(data,1) ~= numel(fits)
    error('Data size must be the same as the number of FITTED objects.');
end
for i = 1:numel(fits)
    fits(i).(name) = data(i,:);
end
end % assign_datafield
