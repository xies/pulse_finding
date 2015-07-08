function H = plot_aligned(cells,name2plot,varargin)
%Plot as a shadedErrorBar the mean specified property of given cellobj
%array. Cells must come from the same embryo (same embryoID).
%
%  USAGE: H = cells.plot_aligned('area')

if numel(unique([cells.embryoID])) > 1
    error('Only cells from the same embryo allowed.')
end

if nargin < 2, name2plot = 'area_sm'; end
if nargin > 2
    color = varargin{1};
else
    color = 'b';
end

time = cells(1).dev_time;
data = cat(2,cells.(name2plot));

H = shadedErrorBar(time,...
    nanmean(data,2), nanstd(data,[],2), ...
    {'color',color},1);

end