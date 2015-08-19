function plot_cells_aligned(pulse,name2plot,varargin)
% Plot average measurement (e.g. area) of a cell across multiple Pulse
% objects (multiple movies).
%
% USAGE: pulse.plot_cells_aligned('area_sm'); -- default behavior

if nargin < 2, name2plot = 'area_sm'; end
if nargin > 2, label2plot = varargin{1}; end

color = hsv(numel(pulse));

for i = 1:numel(pulse)
    
    cellsOI = pulse(i).cells;
    if nargin > 2, cellsOI = cellsOI([cellsOI.label] == label2plot); end
    
    H(i) = cellsOI.plot_aligned(name2plot,color(i,:));
    hold on
    labels{i} = ['Embryo ' num2str(pulse(i).embryoID) ...
        ', dt = ' num2str(pulse(i).input.dt)];
    
end

hold off
ylabel(name2plot)
xlabel('Developmental time')
legend([H.mainLine],labels)

end