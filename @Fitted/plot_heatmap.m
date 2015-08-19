function plot_heatmap(fits,sortname)
%PLOT_HEATMAP Sorts Fitted array by amplitude and then plots the
% concatenated corrected area norm by imagesc.
%
% USAGE: pulse.plot_heatmap % By default sorts by 'amplitude'
%        pulse.plot_heatmap( 'cluster_weight' )
%
% Uses IMAGESC instead of PCOLOR (PCOLOR is werid)

if nargin < 2, sortname = 'amplitude'; end

fits = fits.sort(sortname);

figure

h(1) = subplot(1,4,1:2);
[X,Y] = meshgrid( fits(1).corrected_time, 1:numel(fits));
pcolor( X,Y, cat(1,fits.corrected_myosin) );
%             imagesc( ...
%                 fits(1).corrected_time, ...
%                 1:numel(fits), ...
%                 cat(1,fits.corrected_myosin) );
shading flat; axis tight; colorbar;
title('Myosin intensity')
xlabel('Pulse time (sec)');
axis xy
%             colormap(pmkmp(255))

h(2) = subplot(1,4,3:4);
M = cat(1,fits.corrected_area_norm);
pcolor( X,Y, M);
%             imagesc( ...
%                 fits(1).corrected_time, ...
%                 1:numel(fits), ...
%                 cat(1,fits.corrected_area_norm) );
shading flat; axis tight; colorbar;
caxis( [-8 8] );
title('Local area change');
xlabel('Pulse time (sec)');
axis xy
%             colormap(pmkmp(255))

linkaxes(h)

end % plot_heatmap