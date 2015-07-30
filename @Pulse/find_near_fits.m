function find_near_fits(pulse,neighbor_def)
%FIND_NEAR_FITS Find the number (and fitID) of each fitted
% pulse within an range of time-windows and the first-order
% neighbors, given array of fitted pulses. Results will
% populate its object array.
%
% USAGE: pulse.find_near_fits(neighbor_def)
%
% INPUT: time_windows - 1xN time windows
%        neighbor_def -
%           Fcn handle definition of neighborhood of pulses.
%           Default (ta - tb) > 0.
%           Optionally, can input as a struct with fields:
%               .temporal - temporal window
%               .spatial - spatial window (uses centroid
%               distance instead of graph structure from EDGE)
%
% Note that neighborID returns original EDGE IDs (cellID) and
% not stackIDs.
%
% xies@mit

temp_def = neighbor_def.temporal.def; % definition
time_windows = neighbor_def.temporal.windows;
sp_def = neighbor_def.spatial;

for e = 1:numel(pulse)
    
    % Gather all relevant data into vectors for easy access
    fitsOI = pulse(e).fits; % Gives an array of references to Fitted Of Interest
    cellsOI = pulse(e).cells; % References to CellObj
    nPulse = numel(fitsOI);
    fIDs = [fitsOI.fitID];
    cIDs = [fitsOI.cellID];
    center_frames = [fitsOI.center_frame];
    
    % If center_frame is missing, then populate that field.
    if isempty(center_frames)
        for i = 1:numel(fitsOI)
            I = ...
                findnearest( ...
                fitsOI(i).center,cellsOI(1).dev_time); 
            % H: findnearest Find the nearest numerical 
            % value in an array to a search value 
            % fits.center is a time, (not center of cell)
            if numel(I) > 1
                center_frames(i) = I(1);
            else
                center_frames(i) = I;
            end
        end
    end
    timing = [fitsOI.center];
    
    % Get cellID-cellID spatial conn matrix
    N = cellsOI.get_adjacency_matrix(sp_def);
    
    % Construct fitID-fitID spatial connectivity matrix
    spConn = zeros(nPulse);
    nearCells = zeros(1,nPulse);
    for i = 1:nPulse
        for j = 1:nPulse
            spConn(i,j) = N(cIDs(i),cIDs(j),center_frames(i));
        end
        this_conn = N(cIDs(i),:,center_frames(i));
        nearCells(i) = numel(this_conn(this_conn > 0));
    end
    
    % Get temporal conn matrix based on different timewindow
    % thresholds
    % H: bsxfun(fun,A,B) applies the element-by-element binary operation specified by the function handle fun to arrays A and B
    T = bsxfun(@minus,timing,timing')';
    nearIDs = cell(numel(fitsOI),numel(time_windows)); 
    for i = 1:numel(time_windows)
        
        tempConn = feval(temp_def,T,time_windows(i));   %%%%%%%%%%%%%%%%%   HERE %%%%%
        % H: Multiply two matricies of 1s and 0s
        P = spConn .* tempConn; % full spatiotemporal conn matrix
        
        [I,J] = find(P');
        % Use accumarray with {x} to gather matrix into an cell
        % array
        %H : why is adding a column of ones to J neccessary below ??
        n = accumarray( ...
            cat(2,J,ones(numel(J),1)), fIDs(I), ...
            [size(P,1) 1], @(x) {x} );
        
        % Replace empties with NaNs
        [n{ cellfun(@isempty,n) }] = deal(NaN);
        % Put into larger cell array
        nearIDs(:,i) = n;
        
    end
    
    % Update fitsOI by reference -- no need to return to main workspace
    % Put all data into subsliced Fitted array
    for i = 1:numel(fitsOI)
        fitsOI(i).nearIDs = nearIDs(i,:);
        fitsOI(i).time_windows = time_windows;
        fitsOI(i).neighbor_cells = nearCells(i);
    end
    
end % Loop over all Pulses

end % fin_near_fits