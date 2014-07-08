function MC = reduceMCsimuations(MC_DIR)
% REDUCEMCSIMULATIONS Used to gather together multiple randomization
% tests into a single structure
%
% xies @ mit.edu

% Parse output directory

filelist = what(MC_DIR);
filelist = filelist.mat;

% MC2reduce = cell(1,numel(filelist));

for i = 1:numel(filelist)

	data = load( [MC_DIR '/' filelist{i} ]);
	if isfield(data,'MC');
		MC2reduce(i) = data.MC;
	end

end

MC.random_cell = [MC2reduce(:).random_cell];
MC.empirical = MC2reduce(1).empirical;

MC.neighbor_def = MC2reduce(1).neighbor_def;
MC.time_windows = MC2reduce(1).time_windows;

MC.opt = MC.opt;

end
