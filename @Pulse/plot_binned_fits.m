function plot_binned_fits(pulse)
%Plot error-bar maps of the aligned myosin and aligned area
%change for pulses of different bins

fits = [pulse.fits];
if isempty(fits(1).bin), pulse.bin_fits; end

% Get time vector
x = fits(1).corrected_time;

num_bins = numel(unique(nonans([fits.bin])));
C = pmkmp( num_bins ); % use perceptual map

% iterate through all bin
for i = 1:num_bins
    % Plots either an errorbar or a mean-value (switch
    % depending on plot business)
    fits2bin = fits( [fits.bin] == i);
    
    % myosin subplot
    subplot(2,1,1);
    hold on
    M = nanmean( cat(1, fits2bin.corrected_myosin) );
    %                 shadedErrorBar( x, ...
    %                     nanmean( cat(1, fits2bin.corrected_myosin ) ), ...
    %                     nanstd( cat(1, fits2bin.corrected_myosin ) ), ...
    %                     {'Color',C(i,:)}, 1);
    M = bsxfun(@minus,M,nanmean(M));
    plot(x,M, ...
        'Color',C(i,:),'LineWidth',10);
    xlim([-40 50])
    set(gca,'XTick',[-20 0 20 40]);
    % area subplot
    subplot(2,1,2);
    hold on
    %                 shadedErrorBar( x, ...
    %                     nanmean( cat(1, fits2bin.corrected_area_norm ) ), ...
    %                     nanstd( cat(1, fits2bin.corrected_area_norm ) ), ...
    %                     {'Color',C(i,:)}, 1);
    plot(x,nanmean( cat(1, fits2bin.corrected_area_norm) ), ...
        'Color',C(i,:),'Linewidth',10);
    xlim([-40 50])
    ylim([-4 5])
    set(gca,'XTick',[-20 0 20 40]);
    
end

subplot(2,1,1)
xlabel('Pulse time (sec)')
ylabel('Myosin intensity (a.u.)')
subplot(2,1,2)
xlabel('Pulse time (sec)')
ylabel('\Delta area (\mum^2)')

end % plot_binned_fits