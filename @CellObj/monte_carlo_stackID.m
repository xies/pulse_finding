function [fits,cells] = monte_carlo_stackID(cells,fits,opt)
%MONTE_CARLO_STACKID Randomly exchanges CellObj stackID with
% each other (only fitted + tracked cells). Will also do same
% for FITTED array.
%
% Can specify 'permute' (default) or 'bootstrap'.
%
% USAGE: [fits,cells] =
%               cells.monte_carlo_stackID(fits,'permute')

% Default - do permutation of cells
if nargin < 3, opt = 'permute'; end

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

embryoIDs = unique([fits.embryoID]);
for j = embryoIDs
    
    % extract cells belonging to this embryo (and also has a
    % Fitted object)
    c = cells.get_stackID([fits.get_embryoID(j).stackID]);
    
    % retain original index - for editing
    %                 idx = find( ismember([cells.stackID],[c.stackID]) );
    fidx = { c.fitID };
    % new array of IDs
    sIDs = cat(1,c.stackID);
    cIDs = cat(1,c.cellID);
    
    if strcmpi(opt,'permute')
        % randpermute for permutation
        randIdx = randperm( numel(sIDs) );
    elseif strcmpi(opt,'bootstrap')
        randIdx = randi( numel(sIDs),1,numel(sIDs) );
    else
        error('Unrecognized option - ''permute'' or ''bootstrap''');
    end
    
    % get vector of randomized indices
    sIDs = sIDs( randIdx );
    cIDs = cIDs( randIdx );
    % rewrite cells and fits
    for k = 1:numel(c)
        
        %                     cells( idx(k) ).stackID = sIDs(k);
        %                     cells( idx(k) ).cellID = cIDs(k);
        
        % Rewrite stackID in FITs
        fIDs = [fits(ismember([fits.fitID], fidx{k})).fitID];
        cells( [cells.stackID] == sIDs(k) ).fitID = fIDs;
        
        if any(fIDs)
            fits = fits.set_stackID( fIDs, sIDs(k), cIDs(k) );
        end
    end
end
% sort cells by stackID
[~,I] = sort([cells.stackID],'ascend');
cells = cells(I);

end % bootstrap_stackID