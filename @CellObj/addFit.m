function cellobj = addFit(cellobj,fit)
%@CellObj.addFit Add a fitID to a cell


cellobj.fitID = [cellobj.fitID fit.fitID];
cellobj.num_fits = cellobj.num_fits + 1;
% update gauss curves
cellobj.fit_gausses = ...
    cat(2,cellobj.fit_gausses,...
    lsq_gauss1d([fit.amplitude;fit.center;fit.width],cellobj.fit_time)');
% add new params
cellobj.params = cat(2,cellobj.params,...
    [fit.amplitude;fit.center;fit.width]);
params = cellobj.params;
t = cellobj.fit_time;
% recompute residuals
cellobj.residuals = cellobj.raw ...
    - synthesize_gaussians_withbg(params,t);
P = plot_peak_color(params(:,2:end),t);
cellobj.fit_colorized = P;
end % addFit