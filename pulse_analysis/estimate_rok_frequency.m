%% Auto-correlation analysis
%%
% r_ac - myosin autocorrelation
% rrate_ac - myosin rate autocorrelation

% Fill-in the autoccorrelation matrix (Ncell x Nac):
autocorrelation = mr_ac{1};
frame_rate = 7.5; % sec/frame
x = (0:size(autocorrelation,2)-1)*frame_rate;

%% Find peaks in AC

num_cell = size(autocorrelation,1);
ht = cell(1,num_cell); loc = cell(1,num_cell); freq = cell(1,num_cell);

for j = 1:num_cell
    
    % Find local maxima
    [pk,l] = findpeaks( autocorrelation(j,:),'minpeakdistance',4);
    
    if ~isempty(pk)
        % Gather height (ht) and location (l-1)
        ht{j} = pk; loc{j} = (l-1) * frame_rate;
        % Use intervals in peak location to estimate frequency
        freq{j} = diff([0 loc{j}]);
    end
end

frequency = [freq{:}];
[Nfreq,bins] = hist( frequency, 20);

%% Plot AutoCorr + detected peaks for a single cell

cellID = 40;
plot(x, autocorrelation(cellID,:));
hold on
plot( loc{cellID}, ht{cellID},'ro','MarkerSize',25);
xlabel('Time lag (sec)')
ylabel('Autocorrelation')

%% Plot peak locations for all cells

imagesc(x, 1:num_cell, autocorrelation);
hold on

for i = 1:num_cell
    num_dot = numel(loc{i});
    scatter( loc{i},i*ones(1,num_dot) , 100, 'w','filled')
end
