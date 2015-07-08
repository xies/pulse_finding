%LOAD_EDGE_SCRIPT
%
% Pipeline for loading data from the MATLAB/EDGE-1.06 folder into various file
% structures, including stand-alone arrays (e.g. areas, myosins), as well
% as an array of embryo_structure containing the measurements as fields.
%
% If multiple embryos are loaded, will use the input.tref frame as a
% reference and generate aligned arrays for all data fields (including
% those in embryo_structure)
%
% Important outputs: embryo_stack - 1xN structure containing all loaded
%                      EDGE measurements for N embryos
%
%                    cells - cell-centric structure
%
%                    areas, myosins, etc - T_tot x N_cell num arrays (or
%                      cell-arrays for things like vertices) of the given
%                      data loaded from all embryos
%
%                    input - 1xN structure containing input file/image info
%
%                    num_cells - 1xN array of number of cells loaded from
%                      each embryo
% Feb 2013
% xies@mit.edu

%% Input the filenames and image information

clear input*;

input(1).folder2load = '~/Documents/MATLAB/EDGE-1.06/DATA_GUI/01-30-2012 SqhGap 4/Measurements';
input(1).zslice = 1;
input(1).actual_z = 5;
input(1).tref = 15;
input(1).t0 = 0;
input(1).dt = 6.7;
input(1).um_per_px = .1806;
input(1).X = 1044;
input(1).Y = 400;
input(1).T = 60;
input(1).yref = 32; %um
input(1).embryoID = 1;
input(1).last_segmented = 50;
input(1).fixed = 0;
input(1).ignore_list = [];

input(2).folder2load = '~/Documents/MATLAB/EDGE-1.06/DATA_GUI/01-30-2012 SqhGap 7/Measurements';
input(2).zslice = 2;
input(2).actual_z = 5;
input(2).tref = 30;
input(2).t0 = 0;
input(2).dt = 7.4;
input(2).um_per_px = 0.18921;
input(2).X = 1000;
input(2).Y = 400;
input(2).T = 65; 
input(2).yref = 41; %um
input(2).embryoID = 2;
input(2).last_segmented = 60;
input(2).fixed = 0;
input(2).ignore_list = [];

input(3).folder2load = '~/Documents/MATLAB/EDGE-1.06/DATA_GUI/10-15-2012 SqhGap 1/Measurements';
input(3).zslice = 1;
input(3).actual_z = 4;
input(3).tref = 45;
input(3).t0 = 15;
input(3).ignore_list = [];
input(3).dt = 6.1;
input(3).um_per_px = 0.1732535;
input(3).X = 1000;
input(3).Y = 400;
input(3).T = 80;
input(3).yref = 37; %um
input(3).embryoID = 3;
input(3).last_segmented = 90;
input(3).fixed = 0;
input(3).ignore_list = [];

input(4).folder2load = '~/Documents/MATLAB/EDGE-1.06/DATA_GUI/10-25-2012 SqhGap 1/Measurements';
input(4).zslice = 1;
input(4).actual_z = 4;
input(4).tref = 40;
input(4).t0 = 0;
input(4).ignore_list = [];
input(4).dt = 7.6;
input(4).um_per_px = 0.17969;
input(4).X = 1000;
input(4).Y = 400;
input(4).T = 80;
input(4).yref = 43; %um
input(4).embryoID = 4;
input(4).fixed = 0;
input(4).last_segmented = 65;

input(5).folder2load = '~/Documents/MATLAB/EDGE-1.06/DATA_GUI/11-07-2012 SqhGap 1/Measurements';
input(5).zslice = 2;
input(5).actual_z = 4;
input(5).tref = 90;
input(5).t0 = 0;
input(5).ignore_list = [];
input(5).dt = 7;
input(5).um_per_px = 0.1596724;
input(5).X = 1000;
input(5).Y = 400;
input(5).T = 130;
input(5).yref = 30; %um
input(5).embryoID = 5;
input(5).fixed = 0;
input(5).last_segmented = 125;

% twist

input(6).folder2load = '~/Documents/MATLAB/EDGE-1.06/DATA_GUI/Twist RNAi Series006/Measurements';
input(6).zslice = 1;
input(6).actual_z = 7;
input(6).tref = 1;
input(6).t0 = 0;
input(6).ignore_list = []; %embryo4
input(6).dt = 8;
input(6).um_per_px = 0.141;
input(6).X = 1024; 
input(6).Y = 380;
input(6).T = 100;
input(6).yref = 0; %um
input(6).embryoID = 6;
input(6).fixed = 0;
input(6).last_segmented = 100;

input(7).folder2load = '~/Documents/MATLAB/EDGE-1.06/DATA_GUI/Twist RNAi Series022/Measurements';
input(7).zslice = 1;
input(7).actual_z = 7;
input(7).tref = 10;
input(7).t0 = 0;
input(7).ignore_list = [];
input(7).dt = 8;
input(7).um_per_px = 0.141;
input(7).X = 1024;
input(7).Y = 380;
input(7).T = 70;
input(7).yref = 0; %um
input(7).embryoID = 7;
input(7).fixed = 0;
input(7).last_segmented = 70;

input(8).folder2load = '~/Documents/MATLAB/EDGE-1.06/DATA_GUI/12-04-2013 SqhGap twi RNAi/Measurements';
input(8).zslice = 1;
input(8).actual_z = 6;
input(8).tref = 1;
input(8).t0 = 0;
input(8).ignore_list = [];
input(8).dt = 7.45;
input(8).um_per_px = 0.21255;
input(8).X = 1000;
input(8).Y = 400;
input(8).T = 100;
input(8).fixed = 0;
input(8).yref = 0;
input(8).embryoID = 8;
input(8).last_segmented = 50;

input(9).folder2load = '~/Documents/MATLAB/EDGE-1.06/DATA_GUI/01-15-2014 twi RNAi SqhGap 7/Measurements';
input(9).zslice = 2;
input(9).actual_z = 5;
input(9).tref = 1;
input(9).t0 = 0;
input(9).ignore_list = [];
input(9).dt = 6.961;
input(9).um_per_px = 0.17484;
input(9).X = 1000;
input(9).Y = 400;
input(9).T = 85;
input(9).yref = 0; %um
input(9).embryoID = 9;
input(9).fixed = 0;
input(9).last_segmented = 80;

input(10).folder2load = '~/Documents/MATLAB/EDGE-1.06/DATA_GUI/01-24-2014 twi RNAi SqhGap 5/Measurements';
input(10).zslice = 2;
input(10).actual_z = 7;
input(10).tref = 15;
input(10).t0 = 0;
input(10).ignore_list = [];
input(10).dt = 7.56;
input(10).um_per_px = 0.1844619;
input(10).X = 1000;
input(10).Y = 400;
input(10).T = 100;
input(10).yref = 0; %um
input(10).embryoID = 10;
input(10).fixed = 0;
input(10).last_segmented = 90;

% control injection

input(11).folder2load = '~/Documents/MATLAB/EDGE-1.06/DATA_GUI/02-12-2014 Control SqhGap 1/Measurements';
input(11).zslice = 2;
input(11).actual_z = 7;
input(11).tref = 30;
input(11).t0 = 0;
input(11).ignore_list = [];
input(11).dt = 7.56;
input(11).um_per_px = 0.19385;
input(11).X = 1000;
input(11).Y = 400;
input(11).T = 80;
input(11).yref = 0; %um
input(11).embryoID = 11;
input(11).fixed = 0;
input(11).last_segmented = 70;

input(12).folder2load = '~/Documents/MATLAB/EDGE-1.06/DATA_GUI/02-12-2014 Control SqhGap 4/Measurements';
input(12).zslice = 2;
input(12).actual_z = 7;
input(12).tref = 15;
input(12).t0 = 0;
input(12).ignore_list = [];
input(12).dt = 7.395;
input(12).um_per_px = 0.19385;
input(12).X = 1000;
input(12).Y = 400;
input(12).T = 59;
input(12).yref = 0; %um
input(12).embryoID = 12;
input(12).fixed = 0;
input(12).last_segmented = 59;

input(13).folder2load = '~/Documents/MATLAB/EDGE-1.06/DATA_GUI/02-12-2014 Control SqhGap 2/Measurements';
input(13).zslice = 1;
input(13).actual_z = 8;
input(13).tref = 30;
input(13).t0 = 0;
input(13).ignore_list = [];
input(13).dt = 8;
input(13).um_per_px = 0.141;
input(13).X = 1000;
input(13).Y = 400;
input(13).T = 82;
input(13).yref = 0; %um
input(13).embryoID = 13;
input(13).fixed = 0;
input(13).last_segmented = 80;

input(14).folder2load = '~/Documents/MATLAB/EDGE-1.06/DATA_GUI/Control inject Series001/Measurements';
input(14).zslice = 1;
input(14).actual_z = 8;
input(14).tref = 40;
input(14).t0 = 0;
input(14).ignore_list = [];
input(14).dt = 6.789;
input(14).um_per_px = 0.19385;
input(14).X = 1000;
input(14).Y = 400;
input(14).T = 82;
input(14).yref = 0; %um
input(14).embryoID = 14;
input(14).fixed = 0;
input(14).last_segmented = 80;

input(15).folder2load = '~/Documents/MATLAB/EDGE-1.06/DATA_GUI/Control inject Series002/Measurements';
input(15).zslice = 1;
input(15).actual_z = 8;
input(15).tref = 1;
input(15).t0 = 0;
input(15).ignore_list = [];
input(15).dt = 6.789;
input(15).um_per_px = 0.19385;
input(15).X = 1000;
input(15).Y = 400;
input(15).T = 82;
input(15).yref = 0; %um
input(15).embryoID = 15;
input(15).fixed = 0;
input(15).last_segmented = 50;

% cta

input(16).folder2load = '~/Documents/MATLAB/EDGE-1.06/DATA_GUI/cta_11_10_12_3/Measurements';
input(16).zslice = 1;
input(16).actual_z = 6;
input(16).tref = 16;
input(16).t0 = 0;
input(16).ignore_list = [];
input(16).dt = 6.804;
input(16).um_per_px = 0.2125;
input(16).X = 1000;
input(16).Y = 400;
input(16).T = 100;
input(16).yref = 0; %um
input(16).embryoID = 16;
input(16).fixed = 0;
input(16).last_segmented = 95;
 
% input(17).folder2load = '~/Documents/MATLAB/EDGE-1.06/DATA_GUI/06-03-2014 alphaCat SqhGFP GapmCh 1/Measurements';
% input(17).zslice = 1;
% input(17).actual_z = 6;
% input(17).tref = 15;
% input(17).t0 = 0;
% input(17).ignore_list = [];
% input(17).dt = 7.46;
% input(17).um_per_px = 0.1812;
% input(17).X = 1000;
% input(17).Y = 400;
% input(17).T = 100;
% input(17).yref = 0; %um
% input(17).embryoID = 17;
% input(17).fixed = 0;
% input(17).last_segmented = 100;


%% Rok datasets

clear input_rok
input_rok(1).folder2load = '~/Documents/MATLAB/EDGE-1.06/DATA_GUI/RokGFP_GapCh 140328_3/Measurements';
input_rok(1).zslice = 1;
input_rok(1).actual_z = 5;
input_rok(1).tref = 30;
input_rok(1).t0 = 0;
input_rok(1).ignore_list = []; %embryo4
input_rok(1).dt = 7.5; %embryo4
input_rok(1).um_per_px = 0.2125;
input_rok(1).X = 1000;
input_rok(1).Y = 400;
input_rok(1).T = 103;
input_rok(1).yref = 0; %um
input_rok(1).embryoID = 1;
input_rok(1).fixed = 0;
input_rok(1).last_segmented = 100;

input_rok(2).folder2load = '~/Documents/MATLAB/EDGE-1.06/DATA_GUI/RokGFP_GapCh 140401_2/Measurements';
input_rok(2).zslice = 1;
input_rok(2).actual_z = 5;
input_rok(2).tref = 60;
input_rok(2).t0 = 0;
input_rok(2).ignore_list = []; %embryo4
input_rok(2).dt = 7.5; %embryo4
input_rok(2).um_per_px = 0.2125;
input_rok(2).X = 1000;
input_rok(2).Y = 400;
input_rok(2).T = 120;
input_rok(2).yref = 0; %um
input_rok(2).embryoID = 2;
input_rok(2).fixed = 0;
input_rok(2).last_segmented = 120;

input_rok(3).folder2load = '~/Documents/MATLAB/EDGE-1.06/DATA_GUI/RokGFP GapCherry white RNAi 1/Measurements';
input_rok(3).zslice = 1;
input_rok(3).actual_z = 6;
input_rok(3).tref = 20;
input_rok(3).t0 = 0;
input_rok(3).ignore_list = []; %embryo4
input_rok(3).dt = 9.0877; %embryo4
input_rok(3).um_per_px = 0.141;
input_rok(3).X = 1000; 
input_rok(3).Y = 500;
input_rok(3).T = 30;
input_rok(3).yref = 0; %um
input_rok(3).embryoID = 3;
input_rok(3).fixed = 0;
input_rok(3).last_segmented = 30;

input_rok(4).folder2load = '~/Documents/MATLAB/EDGE-1.06/DATA_GUI/10-15-14 RhoGAP71E HMS E2/Measurements';
input_rok(4).zslice = 1;
input_rok(4).actual_z = 6;
input_rok(4).tref = 27;
input_rok(4).t0 = 0;
input_rok(4).ignore_list = []; %embryo4
input_rok(4).dt = 9; %embryo4
input_rok(4).um_per_px = 0.141;
input_rok(4).X = 1000; 
input_rok(4).Y = 500;
input_rok(4).T = 70;
input_rok(4).yref = 0; %um
input_rok(4).embryoID = 4;
input_rok(4).fixed = 0;
input_rok(4).last_segmented = 30;

% input_gap(2).folder2load = '~/Documents/MATLAB/EDGE-1.06/DATA_GUI/White 10/Measurements';
% input_gap(2).zslice = 1;
% input_gap(2).actual_z = 6;
% input_gap(2).tref = 1;
% input_gap(2).t0 = 0;
% input_gap(2).ignore_list = []; %embryo4
% input_gap(2).dt = 8.3; %embryo4
% input_gap(2).um_per_px = 0.141;
% input_gap(2).X = 1000; 
% input_gap(2).Y = 500;
% input_gap(2).T = 55;
% input_gap(2).yref = 0; %um
% input_gap(2).embryoID = 2;
% input_gap(2).fixed = 0;
% input_gap(2).last_segmented = 55;

msmts2make = {'membranes--basic_2d--area', ...
    'Membranes--vertices--Vertex-y','Membranes--vertices--Vertex-x',...
    'Membranes--basic_2d--Centroid-x','Membranes--basic_2d--Centroid-y',...
    'Membranes--vertices--Identity of neighbors-all', ...
    'Membranes--vertices--Identity of neighbors', ...
    'Myosin--myosin_intensity--Myosin intensity',...
    'Myosin--myosin_intensity--Myosin intensity fuzzy',...
    };

msmts2make_gap = { ...
    'Membranes--basic_2d--area', ...
    'Membranes--basic_2d--Centroid-x','Membranes--basic_2d--Centroid-y',...
    'Myosin--myosin_intensity--Myosin intensity'...
    };

%% Load data (will beep when done)

% EDGEstack = load_edge_data({input.folder2load},msmts2make{:});
% EDGEstack_cta = load_edge_data({input_cta.folder2load},msmts2make{:});
% EDGEstack_princeton = load_edge_data({input_princeton.folder2load},msmts2make_princeton{:});
% EDGEstack_twi = load_edge_data({input_twi.folder2load},msmts2make_twi{:});
EDGEstack_rok = load_edge_data({input_rok.folder2load},msmts2make_gap{:});
beep;

%% Select which stack to load (esp. if there are multiple types of embryos)

in = input_rok;
stack2load = EDGEstack_rok;

%% Load into various stand-alone files (e.g. areas, myosins)

%% membrane

num_embryos = numel(in);

[areas,IDs,dev_time] = extract_msmt_data(stack2load,'area','on',in);
centroids_x = extract_msmt_data(stack2load,'centroid-x','on',in);
centroids_y = extract_msmt_data(stack2load,'centroid-y','on',in);
% neighborID = extract_msmt_data(stack2load,'identity of neighbors-all','off',in);
% vertices_x = extract_msmt_data(stack2load,'vertex-x','off',in);
% vertices_y = extract_msmt_data(stack2load,'vertex-y','off',in);
% majors = extract_msmt_data(stack2load,'major axis','on',input);
% minors = extract_msmt_data(stack2load,'minor axis','on',input);
% % orientations = extract_msmt_data(stack2load,'identity of neighbors','off',input);
% anisotropies = extract_msmt_data(stack2load,'anisotropy','on',input);
% coronal_areas = get_corona_measurement(areas,neighborID);
num_frames = size(areas,1);

areas_sm = smooth2a(areas,1,0);
% rok_sm = smooth2a(squeeze(rok),1,0);
% coronal_areas_sm = smooth2a(coronal_areas,1,0);

areas_rate = -central_diff_multi(areas_sm,dev_time([IDs.which],:));
% anisotropies_rate = central_diff_multi(anisotropies);
% coronal_areas_rate = -central_diff_multi(coronal_areas_sm);

num_cells = zeros(1,num_embryos);
for i = 1:num_embryos
    foo = [IDs.which];
    num_cells(i) = numel(foo(foo==i));
end

%% myosin

myosins = extract_msmt_data(stack2load,'myosin intensity','on',in);
% myosin_ring1 = -extract_msmt_data(stack2load,'ring 1','on',in);
% myosin_ring2 = -extract_msmt_data(stack2load,'ring 2','on',in);
% myosin_ring3 = -extract_msmt_data(stack2load,'ring 3','on',in);
% myosin_inside = extract_msmt_data(stack2load,'inside','on',in);
% myosin_dist2border = extract_msmt_data(stack2load,'distance to border','off',in);
% [myosin_dist2border{cellfun(@isempty,myosin_dist2border)}] = deal(NaN);
% myosin_dist2border = cellfun(@nanmean,myosin_dist2border);
% myosin_size = extract_msmt_data(stack2load,'myosin spot size','off',in);
% myosin_number = extract_msmt_data(stack2load,'number of myosin spots','off',in);
% myosin_fraction = extract_msmt_data(stack2load,'fraction of cell area','on',in);
% myosin_coronal_frac = extract_msmt_data(stack2load,'fraction of coronal area','on',in);
% myosin_connection = extract_msmt_data(stack2load,'# cells connected by myosin','on',in);
% myosin_inertia = extract_msmt_data(stack2load,'moment of inertia','on',in);
% myosin_deviation = extract_msmt_data(stack2load,'deviation from centroid','on',in);
% myosin_span_x = extract_msmt_data(stack2load,'span-x','on',in);
% myosin_span_y = extract_msmt_data(stack2load,'span-y','on',in);
% myosin_perc = extract_msmt_data(stack2load,'# cells connected by myosin','on',in);

myosins_sm = smooth2a(squeeze(myosins),2,0);
% myosins_fuzzy_sm = smooth2a(squeeze(myosins_fuzzy),1,0);
myosins_rate = central_diff_multi(myosins_sm,dev_time([IDs.which],:));
% myosins_rate_fuzzy = central_diff_multi(myosins_fuzzy_sm,1:num_frames);
% coronal_myosins_sm = smooth2a(coronal_myosins,1,0);
% coronal_myosins_rate = central_diff_multi(coronal_myosins_sm);

%% rok

% rok = extract_msmt_data(stack2load,'rok_th intensity','on',in);
% rok_inertia = extract_msmt_data(stack2load,'rok_th - radial moment','on',in);
% rok_dev = extract_msmt_data(stack2load,'rok_th - deviation','on',in);
% rok_rate = central_diff_multi(rok_sm,1:num_frames);

%% Load into embryo-centric structure: embryo_stack, cells

embryo_stack = edge2embryo(stack2load,in,num_cells);
% Put extra data into embryo_stack
for i = 1:num_embryos
    embryo_stack(i).myosin_intensity = myosins(:,[IDs.which] == i);
    embryo_stack(i).area_sm = areas_sm(:,[IDs.which] == i);
    embryo_stack(i).myosin_sm = myosins_sm(:,[IDs.which] == i);
end
cells_raw = embryo2cell(embryo_stack);
