function N = get_adjacency_matrix(cells,method)
%GET_ADJACENCY_MATRIX
% Given cells from a single embryo, return the adjacency matrix
% of the cells defined by identity_of_neighbors_all or by
% spatial threshold of inter-centroid distance.
%
% USAGE: N = cells.get_adjacency_matrix; (default = identity
%                                                    mehtod)
%        N = cells.get_adjacency_matrix(threshold)

% Check that all cells are from same embryo
if numel(unique([cells.embryoID])) > 1,
    error('Cells need to be from same embryo');
end

% Establish which method to use
if nargin < 2
    method.def = 'identity';
end

switch method.def
    case 'identity'
        % Use identity-of-neighbors-all
        nConn = cat(2,cells.identity_of_neighbors_all);
        [T,nCell] = size(nConn);
        
        N = zeros(nCell,nCell,T);
        % Loop through all cells - need to do this b/c there is
        % no native sparse ND matrix support
        for t = 1:T
            n = nConn(t,:);
            for i = 1:nCell
                this_neighbID = n{i};
                if ~isnan(this_neighbID)
                    N(i,this_neighbID(this_neighbID > 0),t) = 1;
                end
            end
        end
        
    case 'window'
        % Calculate all pairwise centroid distances
        nCells = numel(cells);
        T = numel(cells(1).dev_time);
        N = zeros(nCells,nCells,T);
        
        cx = cat(2,cells.centroid_x);
        cy = cat(2,cells.centroid_y);
        
        for t = 1:T
            
            D = squareform( pdist(cat(1,cx(t,:),cy(t,:))') );
            D( logical(eye(size(D,1))) ) = Inf;
            
            if isfield(method,'threshold')
                N(:,:,t) = D < method.threshold;
            else
                N(:,:,t) = D;
            end
            
        end
        
end

% Filter out non-curated cells
% non_curatedIDs = setdiff([cells.cellID],[cells.get_curated.cellID]);
% N(non_curatedIDs,:) = 0;
% N(:,non_curatedIDs) = 0;

end