
% Hierarchical clustering and plots
clear; close all; clc

path_scripts = '/Users/borjablanco/Documents/GitHub/Hyper24RSFC_scripts/scripts/BB_scripts';
addpath(genpath(path_scripts));
path_DOTHUB = '/Users/borjablanco/Documents/GitHub/DOT-HUB_toolbox';
addpath(genpath(path_DOTHUB));
origMeshFileName = '/Volumes/CAM_data/neolab/Rob_Mesh/Robmesh_bb.mshs';
load(origMeshFileName, '-mat')

% Load parcel info to assign the correct labels to parcels
load('/Volumes/CAM_data/neoLAB/group/parcel_info.mat')
load('/Volumes/CAM_data/neoLAB/group/gmSens_parcels.mat')

% Load data and separate into HbO and HbR
% Considering all subjects
load('/Volumes/CAM_data/neoLAB/group/groupICA_allsub_ssr.mat')
hbo = groupICA(:, 1:size(groupICA,2)/2);
hbr = groupICA(:, size(groupICA,2)/2+1:end);

load('/Volumes/CAM_data/neoLAB/group/precisionFC_group_data_ssr.mat')
% For each subject
for j = 1:14
    hbo = hbo_group(:,:,j);
    hbr = hbr_group(:,:,j);
      
    nclusters = 10;
    % Compute hierarchical clustering
    figure
    z_hbo = linkage(hbo', 'ward', 'euclidean');
    tree_hbo = dendrogram(z_hbo, 0, 'ColorThreshold', 300);
    clust_hbo = cluster(z_hbo, 'MaxClust', nclusters);
    set(tree_hbo, 'Linewidth', 2)
    
    figure
    z_hbr = linkage(hbr', 'ward', 'euclidean');
    tree_hbr =dendrogram(z_hbr, 0, 'ColorThreshold', 300);
    clust_hbr = cluster(z_hbr, 'MaxClust', nclusters);
    set(tree_hbr, 'Linewidth', 2)
    close all
    % Assign parcels to labels
    parcel_labels = U(U(:,5)~=0);
    parcel_labels(:,2) = linspace(1, length(parcel_labels), length(parcel_labels));
    parcel_labels(:,3) = clust_hbo;
    parcel_labels(:,4) = clust_hbr;
    
    
    % Plot parcels belonging to each cluster in a different color
    hbo_parcel_cl = gmsens_parcels;
    hbr_parcel_cl = gmsens_parcels;
    
    for i = 1:length(parcel_labels)
        tidx = find(gmsens_parcels==parcel_labels(i,1));
        hbo_parcel_cl(tidx) = parcel_labels(i,3);
        hbr_parcel_cl(tidx) = parcel_labels(i,4);
    end
    
    cmap_parc = distinguishable_colors(nclusters+1);
    cmap_parc(1,:) = 0.9; % to make the rest of the brain grey
    
    fig1 = figure; set(fig1,'units', 'normalized', 'outerposition', [0 0 0.8 0.8], 'Color', [1 1 1]);
    
    subplot(221)
    bb_DOTHUB_plotSurfaceImage (gmSurfaceMesh, hbo_parcel_cl, [180,30], [], cmap_parc)
    caxis([0 nclusters+1]);colorbar off
    
    subplot(222)
    bb_DOTHUB_plotSurfaceImage (gmSurfaceMesh, hbo_parcel_cl, [0,0], [], cmap_parc)
    caxis([0 nclusters+1]);
    colorbar off
    
    subplot(223)
    bb_DOTHUB_plotSurfaceImage (gmSurfaceMesh, hbr_parcel_cl, [180,30], [], cmap_parc)
    caxis([0 nclusters+1]);colorbar off
    
    subplot(224)
    bb_DOTHUB_plotSurfaceImage (gmSurfaceMesh, hbr_parcel_cl, [0,0], [], cmap_parc)
    caxis([0 nclusters+1]);colorbar off
    
    pause
    close all
end  
    








