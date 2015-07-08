function plot_cells_aligned(pulse,name2plot)
% Plot average measurement (e.g. area) of a cell across multiple Pulse
% objects (multiple movies).
%
% USAGE: pulse.plot_cells_aligned('area_sm'); -- default behavior

if nargin < 2, name2plot = 'area_sm'; end

color = hsv(numel(pulse));

for i = 1:numel(pulse)
    
    H(i) = pulse(i).cells.plot_aligned(name2plot,color(i,:));
    hold on
    labels{i} = ['Embryo ' num2str(pulse(i).embryoID) ...
        ', dt = ' num2str(pulse(i).input.dt)];
    
end

hold off
ylabel(name2plot)
xlabel('Developmental time')
legend([H.mainLine],labels)

end