function MC = filter_mc(MC,embryoID)

fieldnames = fields(MC.empirical);
fieldnames= setdiff(fieldnames,{'neighbor_windows'});

for i = 1:numel(fieldnames)
    
    I = ismember(MC.empirical.embryoID,embryoID);
    if isvector(MC.empirical.(fieldnames{i}))
        MC.empirical.(fieldnames{i}) = ensure_column(MC.empirical.(fieldnames{i}));
    end
        
    MC.empirical.(fieldnames{i}) = ...
        MC.empirical.(fieldnames{i})(I,:);
end

N = numel(MC.random_cell);
for n = 1:N
    
    I = ismember(MC.random_cell(n).embryoID,embryoID);
    for i = 1:numel(fieldnames)
        if isvector(MC.random_cell(n).(fieldnames{i}))
            MC.random_cell(n).(fieldnames{i}) = ensure_column(MC.random_cell(n).(fieldnames{i}));
        end
        MC.random_cell(n).(fieldnames{i}) = ...
            MC.random_cell(n).(fieldnames{i})(I,:);
    end
    
end

end