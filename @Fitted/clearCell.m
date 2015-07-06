function fit_array = clearCell(fit_array)
% Clear all records of fits and tracks from a cell
for i = 1:numel(fit_array)
    
    fitOI = fit_array(i);
    % reset all attributes relating to cell
    fitOI.category = [];
    fitOI.bootstrapped = 1;
    fitOI.cellID = [];
    fitOI.stackID = [];
    
    fit_array(i) = fitOI;
    
end
end