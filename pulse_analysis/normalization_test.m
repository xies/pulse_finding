%%

A = zeros(10);

A(1,1) = 3; A(end,1) = 3; A(1,end) = 3; A(end,end) = 3;
A(2:end-1,1) = 5; A(1,2:end-1) = 5; A(2:end-1,end) = 5; A(end,2:end-1) = 5;
A(2:end-1,2:end-1) = 8;

bins = linspace(0,1,5);

%%

Npulse = 15;
Niter = 1000;
N = cell(Niter,8);
M = cell(Niter,8);

for i = 1:Niter
    
    I = randi(89,1,Npulse);
    [X,Y] = find(neighbor_adj(I,:));
    
    nearby = accumarray(X,Y,[numel(I),1],@(x) {x});

    n = cellfun(@(x) numel(x(ismember(x,I))), nearby);
    
    for j = 1:8
        if any(A(I) == j)
            N{i,j} = n(A(I) == j)';
        end
    end
    
%     M{i,1} = flat([n(1,1),n(end,1),n(1,end),n(end,end)]);
%     M{i,2} = [n(2:end-1,1)',n(1,2:end-1),n(2:end-1,end)',n(end,2:end-1)]';
%     M{i,3} = flat(n(2:end-1,2:end-1));
%     
%     N{i,1} = flat([ratio(1,1),ratio(end,1),ratio(1,end),ratio(end,end)]);
%     N{i,2} = [ratio(2:end-1,1)',ratio(1,2:end-1),ratio(2:end-1,end)',ratio(end,2:end-1)]';
%     N{i,3} = flat(ratio(2:end-1,2:end-1));
    
end

%%

[X,G] = make_boxplot_args(cat(1,N{:,1}),cat(1,N{:,2}),cat(1,N{:,3}),{'1','2','3'});
boxplot(X,G)

%%

for i = 1:3
    subplot(3,1,i)
    hist(cat(1,N{:,i}),6);
    xlim([0 3])
end
