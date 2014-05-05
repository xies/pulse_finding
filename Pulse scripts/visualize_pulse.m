%%Visualize_pulse

%% Make movies of individual pulses

% F = make_pulse_movie(fits(101),input,vertices_x,vertices_y,dev_time);
F = fits.movie(101,embryo_stack,cells);
play_movie(F)

%% Make a movie of hand-tracked pulse

G = make_pulse_movie(tracked_pulse(196),input,vertices_x,vertices_y,master_time);

%%
% save movie (to appropriate folder)
% if strcmpi(in(1).folder2load,input_twist(1).folder2load)
%     if ids(cellID).which == 1, var_name = '006'; else var_name = '022'; end
%     movie2avi(f,['~/desktop/edge processed/twist ' var_name '/pulse_movies/pulse_' num2str(pulseID)]);
% elseif strcmpi(in(1).folder2load,input(1).folder2load)
%     if ids(cellID).which == 1, var_name = '4'; else var_name = '7'; end
%     movie2avi(f,['~/desktop/edge processed/embryo '  var_name '/pulse_movies/pulse_' num2str(pulseID)]);
% end

%% Plot individual pulses

ID = 180;

fitsOI = fits_wt.get_cluster(3);
fitsOI = fitsOI.sort('cluster_weight');

t = fits(1).corrected_time;

figure(100),
[ax,h] = plotyy(t,fitsOI(ID).corrected_area_norm,t,fitsOI(ID).corrected_myosin);
set(h,'Color','m')
set(ax(1),'Ycolor','m');
ylabel(ax(1),'Local area response (\mum^2)');
set(ax(2),'Ycolor','g');
ylabel(ax(2),'Myosin intensity (a.u.)');
xlabel('Pulse time (sec)')
fitsOI(ID).fitID
