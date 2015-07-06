function bin_fits(fits,range)
%BIN_FITS Bin fits according to their amplitudes. Quartile binning.

[fits.bin] = deal(NaN);

embryoIDs = unique( [fits.embryoID] );

if nargin < 2, range = 10:10:90; end

for i = embryoIDs
    
    % Get fits in this embryo that are not manually curated
    %                 fits_this_embryo = ...
    %                     fits( [fits.embryoID] == i & ~[fits.manually_added] );
    fits_this_embryo = fits( [fits.embryoID] == i);
    fitIDs_this_embryo = [ fits_this_embryo.fitID ];
    
    % Sort by amplitude
    amps = [ fits_this_embryo.amplitude ];
    [ sorted_amps, sortID ] = sort( amps, 2, 'ascend' );
    
    % Get percentile cutoffs
    cutoffs = prctile( sorted_amps, range );
    cutoffs = [0, cutoffs, max(amps) - eps];
    num_bins = numel(cutoffs) - 1;
    sortedID = fitIDs_this_embryo( sortID );
    
    % Get fitIDs
    binned = cell(1,num_bins);
    for j = 1:num_bins
        if j == numel(cutoffs) - 1
            binned{j} = sortedID( ...
                find(sorted_amps > cutoffs(j),1) : end );
            
        else
            binned{j} = sortedID( ...
                find(sorted_amps > cutoffs(j),1) : ...
                find(sorted_amps > cutoffs(j+1),1) - 1 );
        end
        
    end
    
    % consistency tests
    if numel(unique([binned{:}])) ~= numel(fits_this_embryo) || ...
            ~isempty( setxor([binned{:}],sortedID) )
        display('Something is wrong');
        keyboard;
    end
    
    % Assign bins
    for j = 1:num_bins
        [fits( ismember([fits.fitID],binned{j} ) ).bin] = deal(j);
    end
    
end

end %bin_fits