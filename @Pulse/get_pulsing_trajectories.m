function [adj,nodes] = get_pulsing_trajectories(pulse,varargin)
% GET_PULSING_TRAJECTORIES Construct a graph showing the
% trajectory of a cell through pulse cluster-identity space
% USAGE: [adj,nodes] =
%        cells.get_pulsing_trajectories(fits);
%
% To visualize: wgPlot(adj,nodes)
% Updated: xies@mit Jan 2014

fits = [pulse.fits];

if nargin < 2
    cells = [pulse.cells];
else
    cells = varargin{1};
end

max_fits = nanmax( [cells.num_fits] );
num_clusters = numel(unique([fits.cluster_label]));
adj = zeros(num_clusters*max_fits + 2);

for j = 1:numel(cells)
    
    % Grab current cell, its fits, and fit labels
    this_cell = cells(j);
    this_fits = fits.get_fitID(this_cell.fitID);
    % sort fits by timing
    [~,I] = sort([this_fits.center]);
    states = [this_fits(I).cluster_label];
    
    if ~isempty(states)
        % Origin to first state
        adj( 1,states(1) + 1 ) = adj( 1,states(1) + 1 ) + 1;
        
        for i = 1:numel( states ) - 1
            %States are separated mod(num_clusters)
            beg_index = num_clusters * (i - 1) + states(i) + 1;
            end_index = num_clusters * i + states(i + 1) + 1;
            adj( beg_index, end_index ) = ...
                adj( beg_index, end_index ) + 1;
            
        end
        
        % Last state to end
        adj( end_index, num_clusters*max_fits + 2) = ...
            adj( end_index, num_clusters*max_fits + 2) + 1;
    end
end

x = zeros(num_clusters*max_fits + 2,1);
y = zeros(num_clusters*max_fits + 2,1);

% Construct the nodes
% origin
x(1) = 0; y(1) = (num_clusters + 1)/2;
% ending
x(end) = max_fits + 1; y(end) = (num_clusters + 1)/2;

for i = 1:max_fits
    
    x( (i-1)*num_clusters + 2 : num_clusters*i + 1) = i;
    y( (i-1)*num_clusters + 2 : num_clusters*i + 1) = 1:num_clusters;
    
end
%             k = 1:num_clusters*max_fits;
%             adj = adj(k,k);
nodes = [ x y ];

end % get_pulsing_trajectories