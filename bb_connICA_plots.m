
% NEED TO REORDER MATRICES!!!

% Program to plot in 3D the channels correlations
clear; close all; clc

load('/Volumes/CAM_data/neoLAB/group/parcel_info.mat')
load('/Users/borjablanco/Documents/GitHub/DOT-HUB_toolbox/Utilities/greyJet.mat')
origMeshFileName = '/Volumes/CAM_data/neolab/Rob_Mesh/Robmesh_bb.mshs';
load(origMeshFileName, '-mat')
% Load parcel info to assign the correct labels to parcels
load('/Volumes/CAM_data/neoLAB/group/parcel_info.mat')
load('/Volumes/CAM_data/neoLAB/group/gmSens_mask.mat')
load('/Volumes/CAM_data/neoLAB/group/gmSens_parcels.mat')

load('/Volumes/CAM_data/neoLAB/group/connICA_results_ssr.mat')

edges = size(connICA.IC_comp,1)/2;

hbo_sm = connICA.IC_comp(1:edges,:)';
hbr_sm = connICA.IC_comp(edges+1:2*edges,:)';

hbo_sm = zscore(hbo_sm');
hbr_sm = zscore(hbr_sm');

nparcels = length(unique(gmsens_parcels))-1; % remove parcel of zeros

mask = triu(true(nparcels, nparcels),1);

for nIC=1:39%size(connICA.IC_comp,2)
    
    hbo_vect = hbo_sm(:, nIC);
    % Keep max abs 10%
    max_idx = sort(abs(hbo_vect), 'descend');
    max_val = max_idx(round(length(max_idx)*0.1));
    
    hbo_mat = zeros(nparcels, nparcels);
    hbo_mat(mask) = hbo_vect.*(abs(hbo_vect)>=max_val);
    hbo_mat = hbo_mat + hbo_mat';
    
    hbr_vect = hbr_sm(:, nIC);
    % Keep max abs 10%
    max_idx = sort(abs(hbr_vect), 'descend');
    max_val = max_idx(round(length(max_idx)*0.1));
    hbr_mat = zeros(nparcels, nparcels);
    hbr_mat(mask) = hbr_vect.*(abs(hbr_vect)>=max_val);
    hbr_mat = hbr_mat + hbr_mat';
    
    fig1 = figure; set(fig1, 'units', 'normalized', 'outerposition', [0 0 1 1], 'Color', [1 1 1]);
    subplot(221)
    imagesc(hbo_mat, [-max(abs(hbo_mat(:))) max(abs(hbo_mat(:)))]); axis square
    title('HbO')
    hold on
    % HbO
    % Visual C
    line([0.5 0.5],[0.5 22.5], 'Color', 'k','LineWidth',1.5 )
    line([0.5 22.5],[0.5 0.5], 'Color', 'k','LineWidth',1.5 )
    line([22.5 22.5],[0.5 22.5], 'Color', 'k','LineWidth',1.5 )
    line([0.5 22.5],[22.5 22.5], 'Color', 'k','LineWidth',1.5 )
    
    % Visual per
    line([22.5 30.5],[22.5 22.5], 'Color', 'k','LineWidth',1.5 )
    line([22.5 22.5],[22.5 30.5], 'Color', 'k','LineWidth',1.5 )
    line([30.5 30.5],[22.5 30.5], 'Color', 'k','LineWidth',1.5 )
    line([22.5 30.5],[30.5 30.5], 'Color', 'k','LineWidth',1.5 )
    
    % Dorsal
    line([30.5 43.5],[30.5 30.5], 'Color', 'k','LineWidth',1.5 )
    line([30.5 30.5],[30.5 43.5], 'Color', 'k','LineWidth',1.5 )
    line([43.5 43.5],[30.5 43.5], 'Color', 'k','LineWidth',1.5 )
    line([30.5 43.5],[43.5 43.5], 'Color', 'k','LineWidth',1.5 )
    
    % Salience
    line([43.5 49.5],[43.5 43.5], 'Color', 'k','LineWidth',1.5 )
    line([43.5 43.5],[43.5 49.5], 'Color', 'k','LineWidth',1.5 )
    line([49.5 49.5],[43.5 49.5], 'Color', 'k','LineWidth',1.5 )
    line([43.5 49.5],[49.5 49.5], 'Color', 'k','LineWidth',1.5 )
    
    % Control
    line([49.5 49.5],[49.5 63.5], 'Color', 'k','LineWidth',1.5 )
    line([49.5 63.5],[49.5 49.5], 'Color', 'k','LineWidth',1.5 )
    line([63.5 63.5],[49.5 63.5], 'Color', 'k','LineWidth',1.5 )
    line([49.5 63.5],[63.5 63.5], 'Color', 'k','LineWidth',1.5 )
    
    % Default
    line([63.5 80.5],[63.5 63.5], 'Color', 'k','LineWidth',1.5 )
    line([63.5 63.5],[63.5 80.5], 'Color', 'k','LineWidth',1.5 )
    line([80.5 80.5],[63.5 80.5], 'Color', 'k','LineWidth',1.5 )
    line([63.5 80.5],[80.5 80.5], 'Color', 'k','LineWidth',1.5 )
    
    % Visual C right
    line([80.5 102.5],[80.5 80.5], 'Color', 'k','LineWidth',1.5 )
    line([80.5 80.5],[80.5 102.5], 'Color', 'k','LineWidth',1.5 )
    line([80.5 102.5],[102.5 102.5], 'Color', 'k','LineWidth',1.5 )
    line([102.5 102.5],[80.5 102.5], 'Color', 'k','LineWidth',1.5 )
    
    % Visual per right
    line([102.5 107.5],[102.5 102.5], 'Color', 'k','LineWidth',1.5 )
    line([102.5 102.5],[102.5 107.5], 'Color', 'k','LineWidth',1.5 )
    line([107.5 107.5],[102.5 107.5], 'Color', 'k','LineWidth',1.5 )
    line([102.5 107.5],[107.5 107.5], 'Color', 'k','LineWidth',1.5 )
    
    % Dorsal right
    line([107.5 120.5],[107.5 107.5], 'Color', 'k','LineWidth',1.5 )
    line([107.5 107.5],[107.5 120.5], 'Color', 'k','LineWidth',1.5 )
    line([120.5 120.5],[107.5 120.5], 'Color', 'k','LineWidth',1.5 )
    line([107.5 120.5],[120.5 120.5], 'Color', 'k','LineWidth',1.5 )
    
    % Salience right
    line([120.5 128.5],[120.5 120.5], 'Color', 'k','LineWidth',1.5 )
    line([120.5 120.5],[120.5 128.5], 'Color', 'k','LineWidth',1.5 )
    line([128.5 128.5],[120.5 128.5], 'Color', 'k','LineWidth',1.5 )
    line([120.5 128.5],[128.5 128.5], 'Color', 'k','LineWidth',1.5 )
    
    % Control right
    line([128.5 143.5],[128.5 128.5], 'Color', 'k','LineWidth',1.5 )
    line([128.5 128.5],[128.5 143.5], 'Color', 'k','LineWidth',1.5 )
    line([143.5 143.5],[128.5 143.5], 'Color', 'k','LineWidth',1.5 )
    line([128.5 143.5],[143.5 143.5], 'Color', 'k','LineWidth',1.5 )
    
    % Default right
    line([143.5 154.5],[143.5 143.5], 'Color', 'k','LineWidth',1.5 )
    line([143.5 143.5],[143.5 154.5], 'Color', 'k','LineWidth',1.5 )
    line([154.5 154.5],[143.5 154.5], 'Color', 'k','LineWidth',1.5 )
    line([143.5 154.5],[154.5 154.5], 'Color', 'k','LineWidth',1.5 )
    
    
    
    subplot(222)
    imagesc(hbr_mat, [-max(abs(hbr_mat(:))) max(abs(hbr_mat(:)))]); axis square
    title('HbR')
    hold on
    % HbR
    % Visual C
    line([0.5 0.5],[0.5 22.5], 'Color', 'k','LineWidth',1.5 )
    line([0.5 22.5],[0.5 0.5], 'Color', 'k','LineWidth',1.5 )
    line([22.5 22.5],[0.5 22.5], 'Color', 'k','LineWidth',1.5 )
    line([0.5 22.5],[22.5 22.5], 'Color', 'k','LineWidth',1.5 )
    
    % Visual per
    line([22.5 30.5],[22.5 22.5], 'Color', 'k','LineWidth',1.5 )
    line([22.5 22.5],[22.5 30.5], 'Color', 'k','LineWidth',1.5 )
    line([30.5 30.5],[22.5 30.5], 'Color', 'k','LineWidth',1.5 )
    line([22.5 30.5],[30.5 30.5], 'Color', 'k','LineWidth',1.5 )
    
    % Dorsal
    line([30.5 43.5],[30.5 30.5], 'Color', 'k','LineWidth',1.5 )
    line([30.5 30.5],[30.5 43.5], 'Color', 'k','LineWidth',1.5 )
    line([43.5 43.5],[30.5 43.5], 'Color', 'k','LineWidth',1.5 )
    line([30.5 43.5],[43.5 43.5], 'Color', 'k','LineWidth',1.5 )
    
    % Salience
    line([43.5 49.5],[43.5 43.5], 'Color', 'k','LineWidth',1.5 )
    line([43.5 43.5],[43.5 49.5], 'Color', 'k','LineWidth',1.5 )
    line([49.5 49.5],[43.5 49.5], 'Color', 'k','LineWidth',1.5 )
    line([43.5 49.5],[49.5 49.5], 'Color', 'k','LineWidth',1.5 )
    
    % Control
    line([49.5 49.5],[49.5 63.5], 'Color', 'k','LineWidth',1.5 )
    line([49.5 63.5],[49.5 49.5], 'Color', 'k','LineWidth',1.5 )
    line([63.5 63.5],[49.5 63.5], 'Color', 'k','LineWidth',1.5 )
    line([49.5 63.5],[63.5 63.5], 'Color', 'k','LineWidth',1.5 )
    
    % Default
    line([63.5 80.5],[63.5 63.5], 'Color', 'k','LineWidth',1.5 )
    line([63.5 63.5],[63.5 80.5], 'Color', 'k','LineWidth',1.5 )
    line([80.5 80.5],[63.5 80.5], 'Color', 'k','LineWidth',1.5 )
    line([63.5 80.5],[80.5 80.5], 'Color', 'k','LineWidth',1.5 )
    
    % Visual C right
    line([80.5 102.5],[80.5 80.5], 'Color', 'k','LineWidth',1.5 )
    line([80.5 80.5],[80.5 102.5], 'Color', 'k','LineWidth',1.5 )
    line([80.5 102.5],[102.5 102.5], 'Color', 'k','LineWidth',1.5 )
    line([102.5 102.5],[80.5 102.5], 'Color', 'k','LineWidth',1.5 )
    
    % Visual per right
    line([102.5 107.5],[102.5 102.5], 'Color', 'k','LineWidth',1.5 )
    line([102.5 102.5],[102.5 107.5], 'Color', 'k','LineWidth',1.5 )
    line([107.5 107.5],[102.5 107.5], 'Color', 'k','LineWidth',1.5 )
    line([102.5 107.5],[107.5 107.5], 'Color', 'k','LineWidth',1.5 )
    
    % Dorsal right
    line([107.5 120.5],[107.5 107.5], 'Color', 'k','LineWidth',1.5 )
    line([107.5 107.5],[107.5 120.5], 'Color', 'k','LineWidth',1.5 )
    line([120.5 120.5],[107.5 120.5], 'Color', 'k','LineWidth',1.5 )
    line([107.5 120.5],[120.5 120.5], 'Color', 'k','LineWidth',1.5 )
    
    % Salience right
    line([120.5 128.5],[120.5 120.5], 'Color', 'k','LineWidth',1.5 )
    line([120.5 120.5],[120.5 128.5], 'Color', 'k','LineWidth',1.5 )
    line([128.5 128.5],[120.5 128.5], 'Color', 'k','LineWidth',1.5 )
    line([120.5 128.5],[128.5 128.5], 'Color', 'k','LineWidth',1.5 )
    
    % Control right
    line([128.5 143.5],[128.5 128.5], 'Color', 'k','LineWidth',1.5 )
    line([128.5 128.5],[128.5 143.5], 'Color', 'k','LineWidth',1.5 )
    line([143.5 143.5],[128.5 143.5], 'Color', 'k','LineWidth',1.5 )
    line([128.5 143.5],[143.5 143.5], 'Color', 'k','LineWidth',1.5 )
    
    % Default right
    line([143.5 154.5],[143.5 143.5], 'Color', 'k','LineWidth',1.5 )
    line([143.5 143.5],[143.5 154.5], 'Color', 'k','LineWidth',1.5 )
    line([154.5 154.5],[143.5 154.5], 'Color', 'k','LineWidth',1.5 )
    line([143.5 154.5],[154.5 154.5], 'Color', 'k','LineWidth',1.5 )
    
    
    
    sgtitle(['Component ' num2str(nIC)])
    subplot(223)
    bar(connICA.weights(:,nIC))
    subplot(224)
    bar(connICA.weights(:,nIC))
    colormap jet
    pause
    close all
end



imagesc(corr(connICA.IC_comp), [-1 1]); colormap jet







