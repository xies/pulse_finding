
function cells = update_measurements(cells,embryo_stack)

num_cells = numel(cells);
if numel(cells) ~= size(embryo_stack.area,2)
    error('Number of cells must match')
end

measurements = setdiff(fieldnames(embryo_stack), ...
    {'input','dev_time','dev_frame','num_cell','num_frame',...
    'identity_of_neighbors'});
for i = 1:num_cells
    for m = measurements'
        X = embryo_stack.(m{:});
        cells(i).(m{:}) = X(:,i);
    end
end

end % update_measurements