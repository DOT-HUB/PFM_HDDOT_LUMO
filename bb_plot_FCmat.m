clear; close all; clc
load('/Volumes/CAM_data/neoLAB/group/connICA_allsub_ssr.mat')
load('/Volumes/CAM_data/neoLAB/group/gmSens_parcels.mat')
nparcels = size(groupFC_z,1)/2;
group_hbo = tanh(mean(groupFC_z (1:nparcels, 1:nparcels, :),3));
group_hbr = tanh(mean(groupFC_z (nparcels+1:end, nparcels+1:end, :),3));

% Reorder matrix based on parcels
% Plot parcels according to Schaefer labels based on Yeo parcellation
% Reorder data and FC matrices according to parcel Schaefer
% Load txt with Schaefer labels
fid = fopen('/Volumes/CAM_data/neoLAB/group/Schaefer2018_1000Parcels_17Networks_order.txt');
label_info = textscan(fid,'%f%s%f%f%f%f');
fclose(fid);

% Assign labels to parcels
gm_parc = unique(gmsens_parcels); % find unique parcels, parcel number
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




fig1 = figure; set(fig1, 'units', 'normalized', 'outerposition', [0 0 1 1], 'Color', [1 1 1]);
subplot(1,2,1)
imagesc(group_hbo, [-1 1])
set(gca, 'FontSize', 22)
title ('HbO')
axis square
hold on
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

subplot(1,2,2)
imagesc(group_hbr, [-1 1])
axis square
set(gca, 'FontSize', 22)
title ('HbR')
colormap jet
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

