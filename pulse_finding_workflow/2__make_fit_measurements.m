
%%

clear fit_opts
[fit_opts(1:num_embryos).to_fit] = deal('myosin_intensity');
[fit_opts(1:num_embryos).bg] = deal('on');

[fit_opts(1:num_embryos).left_margin] = deal(6);
[fit_opts(1:num_embryos).right_margin] = deal(10);
[fit_opts(1:num_embryos).nan_thresh] = deal(30);
[fit_opts(1:num_embryos).nan_consec_thresh] = deal(4);
[fit_opts(1:num_embryos).end_tol] = deal(30);

[fit_opts(1:num_embryos).alpha] = deal(0.01);
[fit_opts(1:num_embryos).sigma_lb] = deal(8);
[fit_opts(1:num_embryos).sigma_ub] = deal(35);

% [fit_opts(11).to_fit] = deal('myosin_intensity_fuzzy');

%% Perform fitting

% if ~exist('cells','var'), cells = embryo2cell(embryo_stack); end
[cells_raw,fits_raw] = fit_gaussians(cells_raw,fit_opts);
% save('~/Desktop/fits_raw.mat','fits_raw','cells_raw')

%% Instantiate Pulse

temp_pulses = Pulse(cells,fits) % Sketch

%% sub-set of pulses

% pulse.retrace( fit_opts );
% pulse.measure_fits();
