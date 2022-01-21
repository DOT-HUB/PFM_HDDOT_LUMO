clear; close all; clc
% Load mesh for plots
origMeshFileName = '/Volumes/CAM_data/neolab/Rob_Mesh/Robmesh_bb.mshs';
load(origMeshFileName, '-mat')
path_DOTHUB = '/Users/borjablanco/Documents/GitHub/DOT-HUB_toolbox';
addpath(genpath(path_DOTHUB));

path_scripts = '/Users/borjablanco/Documents/GitHub/Hyper24RSFC_scripts/scripts/BB_scripts';

addpath(genpath(path_scripts));

path_figures = '/Volumes/CAM_data/neoLAB/group/figures';
% Create plots for 60% to 75% according to volcano plots (first 4)
% Load sensitivity - parcel space
load('/Volumes/CAM_data/neoLAB/group/gmSens_parcels.mat');
load('/Volumes/CAM_data/neoLAB/group/parcel_info.mat')
load('/Users/borjablanco/Documents/GitHub/Hyper24RSFC_scripts/scripts/greyJet.mat')
load('/Volumes/CAM_data/neoLAB/group/individual_yeo_networks_mat/yeoNet.mat')

load('/Volumes/CAM_data/neoLAB/group/groupICA_spatial_indivMaps.mat')
load('/Volumes/CAM_data/neoLAB/group/groupICA_spatial_results_ssr.mat')

npca = 3; % up to 70% PCA variance threshold

hbo_pos_all = [];
hbr_pos_all = [];
hbo_neg_all = [];
hbr_neg_all = [];

% Select maps on results' structure
j = sum(groupICA_results.nComp(1:npca-1))+1;
hbo_maps = group_sm_hbo(j:(j+groupICA_results.nComp(npca))-1,:,:);
hbr_maps = group_sm_hbr(j:(j+groupICA_results.nComp(npca))-1,:,:);

% adjust network sign for visualization
%net_sign = [1, -1, 1, 1, 1, 1, 1, -1, -1, 1, 1, ...
%    1, 1, -1, 1, 1, 1]; 65%

net_sign = [1, -1, 1, -1, 1, 1, -1, -1, 1, 1, -1, ...
    1, 1, 1, 1, 1, 1, 1, -1, -1, -1]; % 70%

npar = 154;
nIC = 21;
ica_Tmap_hbo_all = zeros(npar, nIC);
ica_Tmap_hbr_all = zeros(npar, nIC);

cmap_lab = 'vik';
cmap = load(['/Users/borjablanco/Documents/GitHub/Hyper24RSFC_scripts/scripts/BB_scripts/colormats/' cmap_lab]);
cmap = (cmap.vik);

cmap(107:149,:) = repmat(cmap(128,:),43,1);

%cmap = greyJet;

bonf_thresh = 0.05/size(ica_Tmap_hbo_all,1);

% Compute T-maps based on individual groupICA maps from dual regression
for ncomp = 1:groupICA_results.nComp(npca)
    
    gmsens_hbo = zeros(size(gmsens_parcels,1),1);
    ica_map_hbo = squeeze(hbo_maps(ncomp,:,:));
    [h_hbo, ~, ~, stats_hbo] = ttest(ica_map_hbo',0, 'Alpha', bonf_thresh);
    ica_Tmap_hbo = (h_hbo.*stats_hbo.tstat)'*net_sign(ncomp);
    
    for npar = 1:length(ica_Tmap_hbo)
        parc_lab = U(:,5);
        parc_lab = parc_lab(parc_lab~=0);
        yeo_idx = find(gmsens_parcels==parc_lab(npar));
        gmsens_hbo(yeo_idx) = ica_Tmap_hbo(npar);
    end
    
    % Plot HbO
    % Number of subplots
    fig1 = figure; set(fig1, 'units', 'normalized', 'outerposition', [0 0 1 1], 'Color', [1 1 1]);
    subplot(2, 2, 1);
    DOTHUB_plotSurfaceImage(gmSurfaceMesh, gmsens_hbo, [0 0],[], cmap);
    title('HbO')
    colorbar off
    subplot(2, 2, 3);
    DOTHUB_plotSurfaceImage(gmSurfaceMesh, gmsens_hbo, [180 35],[], cmap);
    title(['IC ' num2str(ncomp) ' npca ' num2str(npca)])
    colorbar off

    %HbR
    gmsens_hbr = zeros(size(gmsens_parcels,1),1);
    ica_map_hbr = squeeze(-hbr_maps(ncomp,:,:));
    [h_hbr, p, ~, stats_hbr] = ttest(ica_map_hbr',0, 'Alpha', bonf_thresh);
    ica_Tmap_hbr = (h_hbr.*stats_hbr.tstat)'*net_sign(ncomp);
    
    for npar = 1:length(ica_Tmap_hbr)
        parc_lab = U(:,5);
        parc_lab = parc_lab(parc_lab~=0);
        yeo_idx = find(gmsens_parcels==parc_lab(npar));
        gmsens_hbr(yeo_idx) = ica_Tmap_hbr(npar);
    end
    
    
    ica_Tmap_hbo_all(:,ncomp) = ica_Tmap_hbo;
    ica_Tmap_hbr_all(:,ncomp) = ica_Tmap_hbr;

    % Plot HbR
    subplot(2, 2, 2);
    DOTHUB_plotSurfaceImage(gmSurfaceMesh,gmsens_hbr, [0 0],[], cmap);
    title('HbR')
    colorbar off

    %colorbar('southoutside')
    subplot(2, 2, 4);
    DOTHUB_plotSurfaceImage(gmSurfaceMesh, gmsens_hbr, [180 35],[], cmap);
    title(['IC ' num2str(ncomp) ' npca ' num2str(npca)])
    c = colorbar;
    set(c, 'fontsize', 40, 'Position', [0.88, 0.11, 0.03, 0.5])
    
    pause
    
    cd(path_figures)
    saveas(fig1,['network_' num2str(ncomp) '.tiff'], 'tiff');
    close
    
    
%     % --------------- Change this part by Jaccard idx ------------
%     for i = 1:6; hbo_pos (i) = jaccard(yeoN(:,i), double(gmsens_hbo>0)); end
%     for i = 1:6; hbr_pos (i) = jaccard(yeoN(:,i), double(gmsens_hbr>0)); end
%     for i = 1:6; hbo_neg (i) = jaccard(yeoN(:,i), double(gmsens_hbo<0)); end
%     for i = 1:6; hbr_neg (i) = jaccard(yeoN(:,i), double(gmsens_hbr<0)); end
%     % ------------------------------------------------------------
%      
%     subplot(2,4,5)
%     imagesc(hbo_pos', [0 0.5])
%     title('HbO +')
%     set(gca, 'YtickLabel', {'visual C', 'visual P','control','salience', 'dorsal', 'default'})
%     subplot(2,4,6)
%     imagesc(hbr_pos', [0 0.5])
%     title('HbR +')
%     subplot(2,4,7)
%     imagesc(hbo_neg', [0 0.5])
%     title('HbO -')
%     subplot(2,4,8)
%     imagesc(hbr_neg', [0 0.5])
%     title('HbR -')
%     
%     colormap(greyJet)
%     %pause;
%     
%     hbo_pos_all = [hbo_pos_all; hbo_pos];
%     hbr_pos_all = [hbr_pos_all; hbr_pos ];
%     hbo_neg_all = [hbo_neg_all; hbo_neg ];
%     hbr_neg_all = [hbr_neg_all; hbr_neg ];
end



for i =1:21
    subplot(121)
    imagesc(corr(squeeze(hbo_maps(i,:,:))), [-1 1])
    
    subplot(122)
    imagesc(corr(squeeze(hbr_maps(i,:,:))), [-1 1])
    colormap jet
    pause
    close
end





