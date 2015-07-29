function [fits_bs,cells_bs] = simulate_pulsing(pulse,freqHat)
% Simulates spatially random pulses onto the empirical cell lattice
% using existing FITS as seeds and freqHat to estimate the
% frequency between consecutive pulses within a cell and pcHat
% to estimate the number of pulses within a cell

% set up local stream with /dev/random (UNIX only)
% set up local stream with /dev/random to get random seed
sfd = fopen('/dev/urandom');
seed = fread(sfd, 1, 'uint32');
fclose(sfd);

% check MATLAB version
str = version;
if strcmpi(str(1),'8')
    rng(seed)
else
    stream = RandStream('mt19937ar','Seed',seed); % MATLAB's start-up settings
    RandStream.setDefaultStream(stream);
end

%
fits = [pulse.fits];
cells = [pulse.cells];
fits_bs = fits.copy.clearCell;
cells_bs = cells.copy.clearFitsTracks;

% Repeat for each embryo
for e = 1:numel(pulse)
    
    % sort pulses by their center of timing
    fitsOI = pulse(e).fits.sort('center');
    % get all cells in this embryo THAT HAS PULSES IN EMPIRICAL DATASET
    cellsOI = pulse(e).cells;
    Ncells = numel(cellsOI);
    % clear all pulse data associated with cell
    cellsOI = cellsOI.clearFitsTracks;
    
    % Generate the adjacency matrix as a fcn of time
    A = cellsOI.get_adjacency_matrix;
    
    % keep track with Nframe x Ncell matrix of which cell
    % already pulsed
    already_pulsed = zeros( max([fitsOI.center_frame]),Ncells );
    
    % Loop over all embryos
    for i = 1:numel(fitsOI)
        
        accept = 0;
        this_fit = fitsOI(i);
        frame = this_fit.center_frame;
        % Figure out how many adjacent cells current pulsing cell has
        this_fit.neighbor_cells = sum(A(this_fit.cellID,:,frame));
        
        % TODO: Corner case NaN is center_frame - need to deal with
        % case ... right now just spits out same cell
        if this_fit.neighbor_cells == 0
            accept_move(this_fit,pulse(e).find_cells_with_fit(this_fit));
%             fitsOI(i) = this_pulse;
%             cellsOI( [cellsOI.cellID] == cellOI.cellID) = cellOI;
            already_pulsed(i,cellOI.cellID) = 1;
            continue
        end
        
        % Find the number of neighboing cells to the pulse
        N = cat(2,cellsOI.identity_of_neighbors_all);
        num_neighbors = N(frame,:);
        num_neighbors = cellfun(@(x) numel(x(x > 0)), num_neighbors);
        
        % make sure candidate cells have the same number of
        % neighboring cells (index is the index of cells_in_embryo)
        candidate_range = find(num_neighbors == this_fit.neighbor_cells);
        
        % If there is only a single candidate, automatically
        % accept
        if numel(candidate_range) == 1
            
            cellOI = cellsOI(candidate_range);
            accept_move(this_fit,cellOI);
%             fitsOI(i) = this_pulse;
%             cellsOI(candidate_range) = cellOI;
            already_pulsed(frame,candidate_range) = 1;
            
            if ~isempty([cellsOI.fit_bg]), keyboard; end
            
        else
            
            while ~accept
                
                %                             display(['Randomizing ' num2str(this_pulse.fitID)])
                % Find candidate
                randomID = candidate_range(randi(numel(candidate_range)));
                cellOI = cellsOI(randomID);
                
                % Check that the current cell doesn't already have
                % a pulse at this time
                if already_pulsed(frame,randomID) == 1,
                    accept = 0;
                    if this_fit.fitID == 16028, keyboard; end
                else
                    % additional check for neighbor equality
                    num_neighbors = cellOI.identity_of_neighbors_all{ frame };
                    num_neighbors = numel( num_neighbors( num_neighbors > 0 ) );
                    
                    if num_neighbors ~= this_fit.neighbor_cells
                        keyboard
                        accept = 0;
                    else
                        
                        if cellOI.num_fits == 0
                            % TODO: should we automatically accept if this
                            % is the candidate's first pulse?
                            
                            frame = this_fit.center_frame;
                            
                            % Accept this move
                            % TODO: modify acceptance
                            accept = 1;
                            accept_move(this_fit,cellOI);
%                             fitsOI(i) = this_pulse;
%                             cellsOI(randomID) = cellOI;
                            already_pulsed(frame,randomID) = 1;
                            %                             end
                            
                        else
                            
                            % If there is already a pulse in cell, then
                            % check for interval between pulses
                            interval = this_fit.center - ...
                                max( [fits_bs.get_fitID(cellOI.fitID).center] );
                            
                            % Figure out if input frequency is a histogram or not
                            if isfield(freqHat,'bin') && ~isfield(freqHat,'fun')
                                idx = findnearest(freqHat.bin,interval);
                                p = freqHat.prob(idx);
                            elseif ~isfield(freqHat,'bin') && isfield(freqHat,'fun')
                                p = feval(freqHat.fun,interval);
                            end
                            % How to accept?
                            if rand >= p %rand generates a random uniform number [0,1]
                                % Decline move
                                accept = 0;
                            else
                                % Accept this move
                                accept_move(this_fit,cellOI);
%                                 fitsOI(i) = this_pulse;
%                                 cellsOI(randomID) = cellOI;
                                already_pulsed(frame,randomID) = 1;
                                accept = 1;
                                
                            end % whether to accept based on random number
                            
                        end % Condition on distribution of intervals
                        
                    end % Condition on distribution of number of neighboring cells
                    
                end % make sure cell doesn't already have a pulse
                
            end % while loop for accepting move
            
        end % accept if only 1 candidate
        
        fits_bs( [fits.embryoID] == e ) = fitsOI;
        
    end % Loop over all pulses within embryo
    
    % TODO: Figure out how to re-insert pulse/cell into array
    cells_bs( ismember([cells.stackID],[cellsOI.stackID]) ) ...
        = cellsOI;
    
end % Loop over all embryos

    function accept_move(f,c)
        % Output unnecessary since passed by reference
        c.flag_tracked = 1;
        c.flag_fitted = 1;
%         f.stackID = c.stackID;
        f.cellID = c.cellID;
        f.bootstrapped = 1;
        
        if isnan(c.num_fits), c.num_fits = 0; end
        
        c.num_fits = c.num_fits + 1;
        c.fitID = [c.fitID f.fitID];
    end

end % simulate_pulsing