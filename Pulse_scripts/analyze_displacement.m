%

% Grab appropriate fits + cells
embryoID = 1:5;
% embryoID = 6:10;
c = 3;

fitsOI = fits.get_embryoID(embryoID);
cellsOI = cells.get_embryoID(embryoID);

clear neighbor_defition
neighbor_defition.temporal.def = @(time_diff,tau) (abs(time_diff) < tau);
neighbor_defition.temporal.windows = time_windows;
neighbor_defition.spatial.def = 'identity';

fitsOI = fitsOI.find_near_fits(cellsOI,neighbor_defition);

% Get average displacement of central cells for ALL CELLS
cx = cat(2,cellsOI.centroid_x);
cy = cat(2,cellsOI.centroid_y);
cell_displacements = sqrt( diff(cx).^2 + diff(cy).^2 );
cell_displacements = cat(1,nan(1,size(cell_displacements,2)),cell_displacements);
% Get the angle of displacements
cell_angles = atan2( diff(cy), diff(cx) );
cell_angles = cat(1,nan(1,size(cell_angles,2)),cell_angles);

% Get coronal measurements
coronal_displacements = cellsOI.get_corona_measurement(cell_displacements);
% coronal_displacements = cellfun(@nanmean,coronal_displacements);
coronal_angles = cellsOI.get_corona_measurement(cell_angles);
coronal_angles_difference = cellfun(@minus,coronal_angles,...
    mat2cell(cell_angles,ones(size(cell_angles,1),1), ...
            ones(size(cell_angles,2),1)),'UniformOutput',0);

% % Convert to pulse-centric displacements
self_displacements = fitsOI.get_corrected_measurement(cellsOI,cell_displacements,input);
self_angles = fitsOI.get_corrected_measurement(cellsOI,cell_angles,input);
neighbor_displacements = fitsOI.get_corrected_measurement(cellsOI,coronal_displacements,input);
neighbor_angles = fitsOI.get_corrected_measurement(cellsOI,coronal_angles_difference,input);

% wt(c).R = R(1,2); wt(c).P = P(1,2);
wt(c).vj = [neighbor_displacements{:}];
wt(c).Dtheta = [neighbor_angles{:}];

% twist(c).R = R(1,2); twist(c).P = P(1,2);
% twist(c).vj = [neighbor_displacements{:,7}];
% twist(c).Dtheta = [neighbor_angles{:,7}];

%%

for bin = 1:10
    
    single_bin = fitsOI([fitsOI.bin] == bin);
    
    nearIDs = cat(1,single_bin.nearIDs);
    num_near = cellfun(@(x) numel(x(~isnan(x))), nearIDs);
    
    neighbor_displacements = fitsOI.get_corrected_measurement(cellsOI,coronal_displacements,input);
    neighbor_angles = fitsOI.get_corrected_measurement(cellsOI,coronal_angles_difference,input);
    
    clear coeffs
    for i = 1:numel(single_bin)
        coeffs(i) = nanmean(cos( [neighbor_angles{i,6}]) );
    end
    
    foo = cell(1,16);
    for i = 0:15
        foo{i+1} = coeffs(num_near(:,3) == i);
    end
    
    [foo{cellfun(@isempty,foo)}] = deal(NaN);
    % throw out only 1-element parts
    [foo{cellfun(@(x) numel(x) < 2,foo)}] = deal(NaN);
    
    medians(:,bin) = cellfun(@nanmean,foo);
    
end

medians_wt = medians;
% medians_wt1 = medians;
% medians_wt2 = medians;
% medians_wt3 = medians;
% medians_twist = medians;
% medians_twist1 = medians;
% medians_twist2 = medians;
% medians_twist3 = medians;

%%

subplot(1,2,1),hold on
set(gca,'ColorOrder',pmkmp(10))
plot(0:15,medians_wt),xlim([0 10]),title('WT')
% ylim([-10 15])

subplot(1,2,2),hold all
set(gca,'ColorOrder',pmkmp(10))
plot(0:15,medians_twist),xlim([0 10]),title('twist')
% ylim([-10 15])


%%

figure(1)
subplot(3,1,1),hold on
set(gca,'ColorOrder',pmkmp(10))
plot(0:15,medians_wt1),xlim([0 10]),title('WT')

subplot(3,1,2),hold all
set(gca,'ColorOrder',pmkmp(10))
plot(0:15,medians_wt2),xlim([0 10]),title('twist')

subplot(3,1,3),hold all
set(gca,'ColorOrder',pmkmp(10))
plot(0:15,medians_wt3),xlim([0 10]),title('twist')

xlabel('# neighoring pulses (<60)'),ylabel('Max constriction rate')

%%

figure(3)
subplot(3,1,1),hold on
set(gca,'ColorOrder',pmkmp(10))
plot(0:15,medians_twist1),xlim([0 10]),title('WT')

subplot(3,1,2),hold all
set(gca,'ColorOrder',pmkmp(10))
plot(0:15,medians_twist2),xlim([0 10]),title('twist')

subplot(3,1,3),hold all
set(gca,'ColorOrder',pmkmp(10))
plot(0:15,medians_twist3),xlim([0 10]),title('twist')

xlabel('# neighoring pulses (<60)'),ylabel('Max constriction rate')
