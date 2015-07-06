
function catID = search_catID(pulse,type,pulseID)
%SEARCH_CATID Search for the category_index (catID) of a pulse
% (track/fit). Will use the .category property of the given
% FITTED/TRACK.
%
% USAGE: catID = search_catID(pulse,'fit',fitID)
% 		 catID = search_catID(pulse,'track',trackID)

if strcmpi(type,'fit')
    this = pulse.get_fitID(pulseID);
    ID = 'fitID';
else
    this = pulse.get_trackID(pulseID);
    ID = 'trackID';
end
category = this.category;
% 			curr_cat = pulse.categories.(category);
catID = cellfun(@(x) ismember( pulseID, x ), ...
    {category.(ID)});
catID = find(catID);

end %search_catID