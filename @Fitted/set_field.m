function fits = set_field(fits,fitIDs, fieldname, fieldvalue)
% Find fits in an array with the given fitIDs and set the given
% fieldnames to that fieldvalue
% Will find the dimension of fieldvalue that corresponds to the
% size of the FITTED object array

fitIDs = nonans(fitIDs);

% deal with matrix-valued fieldvalues
if ndims(fieldvalue) > 1
    
    [N,M] = size(fieldvalue);
    which_dim = find( [N,M] == numel(fitIDs) );
    if which_dim == 2, fieldvalue = fieldvalue'; end
    
end

for i = 1:numel(fitIDs)
    hit = [fits.fitID] == fitIDs(i);
    if any(hit)
        fits(hit).(fieldname) = fieldvalue(i,:);
    else
        warning(['FitID ' num2str(fitIDs(i)) ' not found']);
        keyboard
    end
end

end %set_fitID