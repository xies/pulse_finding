function sil = silhouette_width(data,labels,distfun)
%SILHOUETTE_WIDTH Wraper for SILHOUETTE function to get average width.
%
% SYNOPSIS: ws = silhouette_width(data,labels);
%           ws = silhouette_width(data,labels,distfun);
%
% If distfun is not provided, uses euclidiean as default
%
% xies@mit.edu

if nargin == 2, distfun = @nan_eucdist; end

sil = silhouette(data,labels,distfun);
sil = nanmean(sil);

end