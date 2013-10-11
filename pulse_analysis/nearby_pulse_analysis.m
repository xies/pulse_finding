%% Nearby pulse analysis

fitsOI = fits_wt;

%%

time_windows = 10:10:100; % seconds

neighbor_defition = @(central, neighbors, tau) ...
    abs( [neighbors.center] - central.center ) < tau ...
    & ~( neighbors == central ) ...
    & ([neighbors.center] - central.center) >= 0;

fitsOI = fitsOI.find_near_fits(time_windows,neighborID, neighbor_defition);

%%

nearIDs = cat(1,fitsOI.nearIDs);

% Convert to number of pulses
num_near = cellfun(@(x) numel(x(~isnan(x))), nearIDs);

%% MC stackID

entries = {'Ratcheted (stereotyped)','Ratcheted (weak)','Ratcheted (delayed)','Un-ratcheted','Stretched'};
Nboot = 400;

%% Permute stackID

num_member = zeros(1,num_clusters);

random_cell(Nboot).num_near = [];
random_cell(Nboot).origin_labels = [];
random_cell(Nboot).target_labels = [];
random_cell(Nboot).centers = [];

random_pulse(Nboot).num_near = [];
random_pulse(Nboot).origin_labels = [];
random_pulse(Nboot).target_labels = [];
random_pulse(Nboot).centers = [];

for j = 1:Nboot
    
    tic
    % make permutations
    [fits_bs_cell,cells_bs_cell] = cells.bootstrap_stackID(fitsOI);
    % get nearby pulses
    fits_bs_cell = fits_bs_cell.find_near_fits(time_windows,neighborID);
    
    % randomize pulses
    [fits_bs_fit,cells_bs_fit] = fitsOI.bootstrap_stackID(cells);
    fits_bs_fit = fits_bs_fit.find_near_fits(time_windows,neighborID);
    
    % get current cluster for empirical (only once)
    if j == 1
        
        this_nearIDs = cat(1,fitsOI.nearIDs);
        % calculate num-neighbors for each pulse
        num_near = cellfun(@(x) numel(x(~isnan(x))), this_nearIDs);
        % origin labels
        labels = [fitsOI.cluster_label]';
        % get all target labels
        target_labels = ...
            cellfun(@(x) [fitsOI.get_fitID(x).cluster_label], ...
            this_nearIDs,'UniformOutput',0);
        % get centers
        centers = [fitsOI.center];
        
        empirical.num_near = num_near;
        empirical.origin_labels = labels;
        empirical.target_labels = target_labels;
        empirical.centers = centers;
        
    end
    
    % random-cell
    if ~isempty(fits_bs_cell)
        
        nearIDs_cell = cat(1,fits_bs_cell.nearIDs);
        % tabluate num-neighbors for each pulse
        num_near = cellfun(@(x) numel(x(~isnan(x))), nearIDs_cell);
        % get origin labels
        labels = [fits_bs_cell.cluster_label]';
        % break down all target labels
        target_labels = ...
            cellfun(@(x) [fits_bs_cell.get_fitID(x).cluster_label], ...
            nearIDs_cell,'UniformOutput',0);
        centers = [fits_bs_cell.center];
        
        random_cell(j).num_near = num_near;
        random_cell(j).origin_labels = labels;
        random_cell(j).target_labels = target_labels;
        random_cell(j).centers = centers;
        
    end
    
    % random-fit
    if ~isempty(fits_bs_fit)
        
        nearIDs_fit = cat(1,fits_bs_fit.nearIDs);
        % tabluate num-neighbors for each pulse
        num_near = cellfun(@(x) numel(x(~isnan(x))), nearIDs_fit);
        % get origin labels
        labels = [fits_bs_fit.cluster_label]';
        % break down all target labels
        target_labels = ...
            cellfun(@(x) [fits_bs_fit.get_fitID(x).cluster_label], ...
            nearIDs_fit,'UniformOutput',0);
        centers = [fits_bs_fit.center];
        
        random_pulse(j).num_near = num_near;
        random_pulse(j).origin_labels = labels;
        random_pulse(j).target_labels = target_labels;
        random_pulse(j).centers = centers;
    end
    
    T = toc;
    display(['Done with ' num2str(j) ' in ' num2str(T) ' seconds.']);
    
end

MC_twist_pcenter.empirical = empirical;
MC_twist_pcenter.random_cell = random_cell;
MC_twist_pcenter.random_pulse = random_pulse;
MC_twist_pcenter.def = neighbor_defition;

save(['~/Desktop/bootstrap_twist_pcenter_N' num2str(Nboot) '_shorter_sub'], ...
    'MC_twist_pcenter');

%% Select correct timing

% select dataset
MC = MC_wt_mcenter;

window = 3; % neighborhood time window
    
empirical = MC.empirical;
random_cell = MC.random_cell;
random_pulse = MC.random_pulse;

% time bounds
left = [-Inf -Inf 0  60  120 180];
right = [Inf 0    60 120 180 Inf];

for K = 1
    
    num_emp = empirical.num_near;
    num_cell = cat(3,random_cell.num_near);
    num_pulse = cat(3,random_pulse.num_near);
    
    labels_emp = empirical.origin_labels;
    labels_cell = random_cell(1).origin_labels;
    labels_pulse = random_cell(1).origin_labels;
    
    target_emp = empirical.target_labels;
    target_cell = cat(3,random_cell.target_labels);
    target_pulse = cat(3,random_pulse.target_labels);
    
    % filter by time
    filter = @(x) (x.centers > left(K) & x.centers <= right(K));
    
    num_emp = num_emp( filter(empirical), :);
    labels_emp = labels_emp( filter(empirical),: );
    target_emp = target_emp( filter(empirical),:,: );
    
    num_cell = num_cell( filter(random_cell(1)), :, :);
    labels_cell = labels_cell( filter(random_cell(1)) );
    target_cell = target_cell( filter(random_cell(1)), :, :);
    
    num_pulse = num_pulse( filter(random_pulse(1)), :, :);
    labels_pulse = labels_pulse( filter(random_pulse(1)) );
    target_pulse = target_pulse( filter(random_pulse(1)), :, :);
    
    for i = 1:5
    
        % Distribution of means within a behavior
        this_count_cell = squeeze( num_cell( labels_cell == i,window,:) );
        this_count_pulse = squeeze( num_pulse( labels_pulse == i,window,:) );
        this_count_emp = num_emp( labels_emp == i, window);
        
        mean_of_cell = mean(this_count_cell,1);
        mean_of_pulse = mean(this_count_pulse,1);
        mean_of_emp = mean(this_count_emp);
        
        [Nmean_pulse,bins] = hist(mean_of_pulse,30);
        [Nmean_cell,bins] = hist(mean_of_cell,bins);
        
        figure(1)
        subplot(5,1,i);
        h = bar(bins, ...
            cat(1,Nmean_cell/sum(Nmean_cell),Nmean_pulse/sum(Nmean_pulse))', ...
            'LineStyle','None');
        set(h(1),'FaceColor','red'),set(h(2),'FaceColor',[0 100/255 0]);
        vline(mean_of_emp,'b');
        title(entries{i})
        
        if i == 1
            xlabel('Average number of neighbors')
            ylabel('Frequency')
            legend('Random-cell','Random-pulse','Empirical');
        end
        
        zscores_cell(K,i) = ( mean_of_emp - mean(mean_of_cell) ) / std(mean_of_cell);
        zscores_pulse(K,i) = ( mean_of_emp - mean(mean_of_pulse) ) / std(mean_of_pulse);
        
%         breakdown neighbor (tagets)
        for j = 1:6
            
            num_target_emp = cellfun(@(x) numel(x(x == j)), ...
                squeeze( target_emp( labels_emp == i,6,:)) );
            num_target_cell = cellfun(@(x) numel(x(x == j)), ...
                squeeze( target_cell( labels_cell == i,6,:) ));
            num_target_pulse = cellfun(@(x) numel(x(x == j)), ...
                squeeze( target_pulse( labels_pulse == i,6,:) ));
            
            mean_of_emp = mean(num_target_emp);
            mean_of_cell = mean(num_target_cell);
            mean_of_pulse = mean(num_target_pulse);
            
            z_target_cell(j) = (mean_of_emp - mean(mean_of_cell)) ...
                /std(mean_of_cell);
            z_target_pulse(j) = (mean_of_emp - mean(mean_of_pulse)) ...
                /std(mean_of_pulse);
            
        end
%         
%         figure(3)
%         subplot(1,5,(K-1)*5 + i);
%         h = bar(1:6, cat(1,z_target_cell,z_target_pulse)' );
%         set(h(1),'FaceColor','red'),set(h(2),'FaceColor','green');
%         xlim([0 7])
%         ylim([-3 3])
        
    end
    
    figure(2)
    subplot(1,1,K);
    h = bar(1:5,cat(1,zscores_cell(K,:),zscores_pulse(K,:))','LineStyle','none');
    set(h(1),'FaceColor','red'),set(h(2),'FaceColor','green');
    title([ num2str(left(K)) ' < center <= ' num2str(right(K)) ])
    ylim([-2.5 3])
    
end

%% Cross-validation Variance testing

MC = MC_wt_pcenter;

Nboot = numel(MC.random_cell);

window = 6;
Ncv = 1000;
subsamples = [10 50 75 100 250]; % out of Nboots

zscores_cell = zeros(numel(subsamples),100,5);
zscores_pulse = zeros(numel(subsamples),100,5);

for k = 1:numel(subsamples)
    for j = 1:Ncv
        
        IDs = randi(Nboot, 1, subsamples(k));
        random_cell = MC.random_cell( IDs );
        random_pulse = MC.random_pulse( IDs );
        
        num_emp = MC.empirical.num_near;
       
        num_cell = cat(3,random_cell.num_near);
        num_pulse = cat(3,random_pulse.num_near);
        
        labels_emp = MC.empirical.origin_labels;
        labels_cell = random_cell(1).origin_labels;
        labels_pulse = random_pulse(1).origin_labels;
        
        
        for i = 1:num_clusters
            
            % breakdown by cluster behavior
            this_count_emp = num_emp( labels_emp == i, window);
            this_count_cell = squeeze( ...
                num_cell( labels_cell == i, window, :) );
            this_count_pulse = squeeze( ...
                num_pulse( labels_pulse == i, window, :) );
            
            mean_of_cell = mean(this_count_cell,1);
            mean_of_pulse = mean(this_count_pulse,1);
            
            zscores_cell(k,j,i) = ...
                (mean(this_count_emp) - mean(mean_of_cell)) ./ std(mean_of_cell);
            zscores_pulse(k,j,i) = ...
                (mean(this_count_emp) - mean(mean_of_pulse)) ./ std(mean_of_pulse);
            
        end
        
    end
    k
end

Z = cat(4,zscores_cell,zscores_pulse);

for i = 1:5
    subplot(5,1,i)
    errorbar(subsamples,mean(Z(:,:,i,1),2),std(Z(:,:,i,1),[],2),'r-'),hold on
    errorbar(subsamples,mean(Z(:,:,i,2),2),std(Z(:,:,i,2),[],2),'g-')
    title(entries{i});
end
