clear; close all; clc
% Load mesh for plots
origMeshFileName = '/Volumes/CAM_data/neolab/Rob_Mesh/Robmesh_bb.mshs';
load(origMeshFileName, '-mat')
path_DOTHUB = '/Users/borjablanco/Documents/GitHub/DOT-HUB_toolbox';
addpath(genpath(path_DOTHUB));

path_scripts = '/Users/borjablanco/Documents/GitHub/Hyper24RSFC_scripts/scripts/BB_scripts';
addpath(genpath(path_scripts));
% Create plots for 60% to 75% according to volcano plots (first 4)
% Load sensitivity - parcel space
load('/Volumes/CAM_data/neoLAB/group/gmSens_parcels.mat');
load('/Volumes/CAM_data/neoLAB/group/parcel_info.mat')
load('/Users/borjablanco/Documents/GitHub/Hyper24RSFC_scripts/scripts/greyJet.mat')
load('/Volumes/CAM_data/neoLAB/group/individual_yeo_networks_mat/yeoNet.mat')

load('/Volumes/CAM_data/neoLAB/group/groupICA_indivMaps.mat')
load('/Volumes/CAM_data/neoLAB/group/groupICA_results_ssr.mat')

npca = 5; % up to 80% PCA variance threshold

hbo_pos_all = [];
hbr_pos_all = [];
hbo_neg_all = [];
hbr_neg_all = [];

for i = 1:npca
    
    % Select maps on results' structure
    j = sum(groupICA_results.nComp(1:i-1))+1;
    hbo_maps = group_sm_hbo(j:(j+groupICA_results.nComp(i))-1,:,:);
    hbr_maps = group_sm_hbr(j:(j+groupICA_results.nComp(i))-1,:,:);
    
    % Compute T-maps based on individual groupICA maps from dual regression
    for ncomp = 1:groupICA_results.nComp(i)

        gmsens_hbo = zeros(size(gmsens_parcels,1),1);
        ica_map_hbo = squeeze(hbo_maps(ncomp,:,:));
        [h_hbo, ~, ~, stats_hbo] = ttest(ica_map_hbo',0, 'Alpha', 3.2468e-04);
        ica_Tmap_hbo = (h_hbo.*stats_hbo.tstat)';
        
        for npar = 1:length(ica_Tmap_hbo)
            parc_lab = U(:,5);
            parc_lab = parc_lab(parc_lab~=0);
            yeo_idx = find(gmsens_parcels==parc_lab(npar));
            gmsens_hbo(yeo_idx) = ica_Tmap_hbo(npar);
        end
        
        % Plot HbO
        % Number of subplots
        fig1 = figure; set(fig1, 'units', 'normalized', 'outerposition', [0 0 0.5 0.7], 'Color', [1 1 1]);
        subplot(2, 4, 1);
        DOTHUB_plotSurfaceImage(gmSurfaceMesh, gmsens_hbo, [0 0],[], greyJet);
        title('HbO')
        colorbar off
        subplot(2, 4, 2);        
        DOTHUB_plotSurfaceImage(gmSurfaceMesh, gmsens_hbo, [180 35],[], greyJet);
        title(['IC ' num2str(ncomp) ' npca ' num2str(i)])
        colorbar off
        
        %HbR
        gmsens_hbr = zeros(size(gmsens_parcels,1),1);
        ica_map_hbr = squeeze(-hbr_maps(ncomp,:,:));
        [h_hbr, ~, ~, stats_hbr] = ttest(ica_map_hbr',0, 'Alpha', 3.2468e-04);
        ica_Tmap_hbr = (h_hbr.*stats_hbr.tstat)';
        
        for npar = 1:length(ica_Tmap_hbr)
            parc_lab = U(:,5);
            parc_lab = parc_lab(parc_lab~=0);
            yeo_idx = find(gmsens_parcels==parc_lab(npar));
            gmsens_hbr(yeo_idx) = ica_Tmap_hbr(npar);
        end
        
        % Plot HbR
        subplot(2, 4, 3);
        DOTHUB_plotSurfaceImage(gmSurfaceMesh,gmsens_hbr, [0 0],[], greyJet);
        title('HbR')
        colorbar off
        subplot(2, 4, 4);        
        DOTHUB_plotSurfaceImage(gmSurfaceMesh, gmsens_hbr, [180 35],[], greyJet);
        title(['IC ' num2str(ncomp) ' npca ' num2str(i)])
        colorbar off 
        
        hbo_pos = corr(yeoN, gmsens_hbo>0);
        hbr_pos = corr(yeoN, gmsens_hbr>0);
        hbo_neg = corr(yeoN, gmsens_hbo<0);
        hbr_neg = corr(yeoN, gmsens_hbr<0);
      
        subplot(2,4,5)
        imagesc(hbo_pos, [-0.75 0.75])
        title('HbO +')
        set(gca, 'YtickLabel', {'visual C', 'visual P','control','salience', 'dorsal', 'default'})
        subplot(2,4,6)
        imagesc(hbr_pos, [-0.75 0.75])
        title('HbR +')
        subplot(2,4,7)
        imagesc(hbo_neg, [-0.75 0.75])
        title('HbO -')
        subplot(2,4,8)
        imagesc(hbr_neg, [-0.75 0.75])
        title('HbR -')

        colormap(greyJet)
              
        hbo_pos_all = [hbo_pos_all, hbo_pos];
        hbr_pos_all = [hbr_pos_all, hbr_pos ];
        hbo_neg_all = [hbo_neg_all, hbo_neg ];
        hbr_neg_all = [hbr_neg_all, hbr_neg ];
    end
        
    pause
    close all
end

figure
% Find matching with Yeo networks
subplot(411)
imagesc(hbo_pos_all>0.4, [-1 1])
line([14.5 14.5], [0 6.5], 'LineWidth', 1)
line([31.5 31.5], [0 6.5], 'LineWidth', 1)
line([52.5 52.5], [0 6.5], 'LineWidth', 1)
line([77.5 77.5], [0 6.5], 'LineWidth', 1)
colormap jet
set(gca, 'YtickLabel', {'visual C', 'visual P','control','salience', 'dorsal', 'default'})

subplot(412)
imagesc(-hbo_neg_all<-0.4, [-1 1])
line([14.5 14.5], [0 6.5], 'LineWidth', 1)
line([31.5 31.5], [0 6.5], 'LineWidth', 1)
line([52.5 52.5], [0 6.5], 'LineWidth', 1)
line([77.5 77.5], [0 6.5], 'LineWidth', 1)
colormap jet
set(gca, 'YtickLabel', {'visual C', 'visual P','control','salience', 'dorsal', 'default'})

subplot(413)
imagesc(hbr_pos_all>0.4, [-1 1])
line([14.5 14.5], [0 6.5], 'LineWidth', 1)
line([31.5 31.5], [0 6.5], 'LineWidth', 1)
line([52.5 52.5], [0 6.5], 'LineWidth', 1)
line([77.5 77.5], [0 6.5], 'LineWidth', 1)
colormap jet
set(gca, 'YtickLabel', {'visual C', 'visual P','control','salience', 'dorsal', 'default'})

subplot(414)
imagesc(-hbr_neg_all<-0.4, [-1 1])
line([14.5 14.5], [0 6.5], 'LineWidth', 1)
line([31.5 31.5], [0 6.5], 'LineWidth', 1)
line([52.5 52.5], [0 6.5], 'LineWidth', 1)
line([77.5 77.5], [0 6.5], 'LineWidth', 1)
colormap jet
set(gca, 'YtickLabel', {'visual C', 'visual P','control','salience', 'dorsal', 'default'})














