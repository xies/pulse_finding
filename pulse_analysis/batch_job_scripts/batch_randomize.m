%BATCH_RANDOMIZE
% To be used from a master shell script for use in batch jobs.

function batch_simulate_pulsing(INPUT_MAT_FILE,OUT_FILENAME,txtfile)

[DIR_STR,name,ext] = fileparts(INPUT_MAT_FILE);

if ~strcmpi(ext,'.mat'),
	error('Needs a .mat file as input.');
end

disp(['Loading input pulse file: ' INPUT_MAT_FILE]);
load(INPUT_MAT_FILE);

time_windows = 10:10:100; % seconds

neighbor_definition.temporal.def = @(time_diff,tau) (time_diff < tau & time_diff > 0);
neighbor_definition.temporal.windows = time_windows;
neighbor_definition.spatial.def = 'identity';

fits = fits.find_near_fits(cells,neighbor_definition);
o.Nboot = 100;
o.timewindows = time_windows;
o.neighbor_def = neighbor_definition;
o.monte_carlo = 'simulation';
o.filter = 'on';

o.savepath = [OUT_FILENAME]

if strcmpi(txtfile,'on')

    fileID = fopen([OUT_FILENAME 'config.txt'],'w');
    
    % print embryoIDs
    embIDs = unique([fits.embryoID]);
    fprintf(fileID,'%s\t', 'EmbryoIDs');
    fprintf(fileID,'%d ', embIDs);
    fprintf(fileID,'\n');
    
    % print simulation conditions
    fprintf(fileID,'%s\t%d\n', 'NumIter',o.Nboot);
    fprintf(fileID,'%s\t%s\n', 'Spatial neighbor def',neighbor_definition.spatial.def);
    fprintf(fileID,'%s\t%s\n', 'Randomization type',o.monte_carlo);
    
    fclose(fileID);
    
    display(['Log file written to ' OUT_FILENAME 'config.txt']);
    
end

disp('Randomization analysis...')
MC = monte_carlo_pulse_location(fits,cells, o);

end
