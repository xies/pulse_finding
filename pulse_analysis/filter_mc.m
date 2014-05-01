function MC = filter_mc(MC,I)

fieldnames = fields(MC.empirical);
fieldnames= setdiff(fieldnames,{'neighbor_windows','near_angle'});

for i = 1:numel(fieldnames)
    if isvector(MC.empirical.(fieldnames{i}))
        MC.empirical.(fieldnames{i}) = ensure_column(MC.empirical.(fieldnames{i}));
    end
    MC.empirical.(fieldnames{i}) = ...
        MC.empirical.(fieldnames{i})(I,:);
end

N = numel(MC.random_cell);
for n = 1:N
    
    for i = 1:numel(fieldnames)
        if isvector(MC.random_cell(n).(fieldnames{i}))
            MC.random_cell(n).(fieldnames{i}) = ensure_column(MC.random_cell(n).(fieldnames{i}));
        end
        MC.random_cell(n).(fieldnames{i}) = ...
            MC.random_cell(n).(fieldnames{i})(I,:);
    end
    
end
end