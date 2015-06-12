function fits = get_fitID(fit_array,fitID)
% Find the FIT(s) with the given fitID(s)
% USAGE: fitsOI = fits.get_fitID(fitID)
%             fits = fit_array( ismember([ fit_array.fitID ], fitID) );
if iscell(fitID), fits = []; return; end
fitID = nonans(fitID);
fits = Fitted;

if numel(fitID) > 0
    fits(numel(fitID)) = Fitted;
    for i = 1:numel(fitID)
        hit = [fit_array.fitID] == fitID(i);
        if any( hit )
            fits(i) = fit_array( hit );
        end
    end
end

% ger rid of empty ones
fits( cellfun(@isempty,{fits.fitID}) ) = [];
end %get_fitID
