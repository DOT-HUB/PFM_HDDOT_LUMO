clear; close all; clc

path_DOTHUB = '/Users/borjablanco/Documents/GitHub/DOT-HUB_toolbox';
addpath(genpath(path_DOTHUB));

% Load sensitivity - parcel space
load('/Volumes/CAM_data/neoLAB/group/gmSens_parcels.mat');

% Plot parcels according to Schaefer labels based on Yeo parcellation
% Reorder data and FC matrices according to parcel Schaefer
% Load txt with Schaefer labels
fid = fopen('/Volumes/CAM_data/neoLAB/group/Schaefer2018_1000Parcels_17Networks_order.txt');
label_info = textscan(fid,'%f%s%f%f%f%f');
fclose(fid);

% Assign labels to parcels
gm_parc = unique(gmsens_parcels); % find unique parcels
gm_parc(gm_parc==0) = []; % Remove parcel of 0s

% Link parcels to Yeo labels
parc_idx = [];
for i = 1:length(gm_parc)
   parc_idx(i) = find(label_info{:,1}==gm_parc(i));
end

parc_label = label_info{1,2};
parc_label = parc_label(parc_idx);

% Rename to make it shorter
for i = 1:length(parc_label)
    temp_label = parc_label{i,1};
    temp_label = temp_label(1:18);
    parc_label{i,1} = temp_label;
end

[a,b, net_labels] = unique(parc_label);
parc_idx = parc_idx';
parc_idx(:,2) = net_labels;

% Load mesh for plots
origMeshFileName = '/Volumes/CAM_data/neolab/Rob_Mesh/Robmesh_bb.mshs';
load(origMeshFileName, '-mat')

gmsens_yeo = gmsens_parcels;
for i = 1:length(gmsens_yeo)
    atemp = find (gm_parc(:,1) == gmsens_yeo(i));
    if ~isempty(atemp)
        gmsens_yeo(i) = parc_idx(atemp, 2); 
    end
end

% Rename for plots
for i = 1:length(a)
    temp_label = a{i};
    temp_label = temp_label(12:end);
    a{i} = temp_label;
end

aplots = [{' '}; a];

% PLOT PARCELLATION
cmap_parc = distinguishable_colors(13);
cmap_parc(1,:) = 0.9; % to make the rest of the brain grey
fig1 = figure; set(fig1, 'units', 'normalized', 'outerposition', [0 0 1 1], 'Color', [1 1 1]);
subplot(121)
DOTHUB_plotSurfaceImage (gmSurfaceMesh, gmsens_yeo, [0,0], [], cmap_parc)
colorbar off
caxis([0 13])
subplot(122)
DOTHUB_plotSurfaceImage (gmSurfaceMesh, gmsens_yeo, [180,30], [], cmap_parc);
colorbar off
caxis([0 13])
cbh = colorbar ; %Create Colorbar
cbh.Ticks = linspace(0, 13, 14)+0.5 ; %Create 13 ticks from zero to 13
cbh.TickLabels = aplots ;    %Replace the labels of these 8 ticks with the numbers 1 to 8
cbh.TickLabelInterpreter = 'none' ;
cbh.FontSize = 20;


 
 
 
 
 
 
