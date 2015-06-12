function obj = find_nearest_object(pulse,obj_type,cx,cy,ct)
% Use centroids, and mean dev_time to robustly find the
% specified object: 'track' or 'fit'

c = pulse.cells;
cframe = findnearest(ct,c(1).dev_time);
if numel(cframe) > 1, cframe = cframe(1); end
x = cat(2,c.centroid_x);
y = cat(2,c.centroid_y);

x = x(cframe,:); y = y(cframe,:);
d = (x-cx).^2 + (y-cy).^2;

[~,which] = min(d);
stackID = c(which).stackID;

switch obj_type
    case 'track'
        
        obj = pulse.tracks.get_stackID(stackID);
        which = findnearest(cellfun(@nanmean,{obj.dev_time}),ct);
        obj = obj(which);
        
    case 'fit'
        
        obj = pulse.fits.get_stackID(stackID);
        which = findnearest(cellfun(@nanmean,{obj.dev_time}),ct);
        obj = obj(which);
        
    otherwise
        error('Can only take ''type'' or ''fit'' as object type')
end

if isempty(obj), obj = []; end

end