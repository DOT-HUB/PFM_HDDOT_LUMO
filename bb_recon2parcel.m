% ----------------------------
% BB Parcellation script
% 11-02-2021
% ----------------------------
clear; close all; clc

path_scripts = '/Users/borjablanco/Documents/GitHub/Hyper24RSFC_scripts/scripts/BB_scripts';
addpath(genpath(path_scripts));

path_DOTHUB = '/Users/borjablanco/Documents/GitHub/DOT-HUB_toolbox';
addpath(genpath(path_DOTHUB));

path_homer2 = '/Users/borjablanco/Documents/MATLAB/toolboxes/homer2';
addpath(genpath(path_homer2));

path_sd3d = '/Volumes/CAM_data/neolab/SD3D_files';
path_prepro = '/Volumes/CAM_data/neolab/Rest_nirs_files';

% Load mesh and parcellated mesh labels
gmSurfaceMeshLabels = load ('/Volumes/CAM_data/neoLAB/Rob_Mesh/gmSurfaceMeshLabels_Schaefer.mat', '-mat');
gmSurfaceMeshLabels = gmSurfaceMeshLabels.gmSurfaceMeshLabels;

load('/Users/borjablanco/Documents/GitHub/DOT-HUB_toolbox/Utilities/greyJet.mat')

% Create sensitivity mask for each participant
origMeshFileName = '/Volumes/CAM_data/neolab/Rob_Mesh/Robmesh_bb.mshs';
load(origMeshFileName, '-mat')
cd(path_sd3d)
jac_name = dir('*.jac');
rmap_name = dir('*.rmap');
SD_name = dir('*.SD3D');

% % Compute sensitivity profile if not already computed
% GMsens_group = zeros(length(jac_name), length(gmSurfaceMeshLabels));
% GMsens_group_binary = zeros(length(jac_name), length(gmSurfaceMeshLabels));
% path_toast = '/Users/borjablanco/Documents/MATLAB/toolboxes/toast';
% addpath(genpath(path_toast));
% path_toast2 = '/Users/borjablanco/Documents/MATLAB/toolboxes/toastpp-2.0.2';
% addpath(genpath(path_toast2));
% cd(path_prepro)
% prepro_name = dir('*ssr.prepro');
% 
% % Compute GM sensitivity profile to create a group mask
% thresh = 0.05; % why this number? arbitrary, but anything between 0.01 and 0.1
% 
% for nsub = 1:length(prepro_name)
%     
%     % This SD3D file needs to be the updated version after excluding
%     % channels? that's what it says on DOTHUB_MakeGMSensitivityMap
%     cd(path_prepro)
%     load(prepro_name(nsub).name, '-mat');
%     SD = SD3D;
%     cd(path_sd3d)   
%     [J_vol_wav1,J_vol_wav2, J_GM_wav1, J_GM_wav2] = JGMtoJVol(jac_name(nsub).name,rmap_name(nsub).name,origMeshFileName);
% 
%     [GMSensitivity, GMmask, J_GM_wav1_norm, J_GM_wav2_norm] = DOTHUB_MakeGMSensitivityMap(J_GM_wav1,J_GM_wav2,J_vol_wav1, J_vol_wav2,SD,thresh);
%     
%     GMsens_group(nsub,:) = GMSensitivity;
%     GMsens_group_binary(nsub,:) = GMmask;
%     
% end
% 
% % Save GM sens

% Load mask if already computed
load('/Volumes/CAM_data/neoLAB/GMsens_group_05threshold.mat')
load('/Volumes/CAM_data/neoLAB/GMsens_group_binary_05threshold.mat')

% Plot average sensitivity
figure
subplot(121)
DOTHUB_plotSurfaceImage (gmSurfaceMesh, mean(GMsens_group),[0,30],[],greyJet);
colorbar off 
subplot(122)
DOTHUB_plotSurfaceImage (gmSurfaceMesh, mean(GMsens_group),[180,30],[],greyJet);
colorbar off

% Plot sensitivity binary mask - how many subjects 
figure
subplot(121)
DOTHUB_plotSurfaceImage (gmSurfaceMesh, sum(GMsens_group_binary),[0,30],[],jet(28));
colorbar off 
subplot(122)
DOTHUB_plotSurfaceImage (gmSurfaceMesh, sum(GMsens_group_binary),[180,30],[],jet(28));
colorbar off

% PLOT PARCELLATION
cmap_parc = distinguishable_colors(max(gmSurfaceMeshLabels)*4);
figure
subplot(121)
bb_DOTHUB_plotSurfaceImage(gmSurfaceMesh, gmSurfaceMeshLabels, [0,30], [], cmap_parc)
colorbar off
subplot(122)
bb_DOTHUB_plotSurfaceImage (gmSurfaceMesh, gmSurfaceMeshLabels,[180,30], [], cmap_parc);
colorbar off

% Find nodes where the setup has enough sensitivity
% thresh = 0.1; % keep values with X% of the max sensitivity
% max_sens = max(mean(GMsens_group)); % mean across subjects
% sens_idx = (mean(GMsens_group)/max_sens)>=thresh; % create mask
nthresh = 12; % 85%, = min number of sessions sensitive to a region
sens_idx = (sum(GMsens_group_binary)>=nthresh); % all subjects

% Remove cerebellum areas
z_val = 108; % value on the z axis
idx_cbl = find(gmSurfaceMesh.node(:,3)<z_val);
sens_idx(idx_cbl) = 0; % add cerebellum indexes to mask

% Mask parcels with no sensitivity or cerebellum
Jmask = sens_idx.*gmSurfaceMeshLabels';
mask_parcels = unique(Jmask);

% Compute number of nodes inside each parcel after masking
[uv,~,idx] = unique(Jmask);
n = accumarray(idx(:),1);
%n(find(uv==0)) = [];

% Check proportion of nodes included with respect to original value
[U(:,1),~,idx_orig] = unique(gmSurfaceMeshLabels);
U(:,2) = accumarray(idx_orig(:),1);
[r,c] = find(U(:,1)==uv);
U(:,3) = 0;
U(r, 3) = n(2:end); % exclude 0s
% Check proportion
U(:,4) = U(:,3)./U(:,2);

% Mask parcels that retained less than certain % of nodes
node_prop = 0.5;
U(:,5) = (U(:,4)>node_prop).*U(:,1);

% Apply sensitivity to parcellation
% Exclude parcels not included in Jmask
gmsens_parcels = zeros(length(gmSurfaceMeshLabels),1);

for i = 1:length(U)
    idx = find(Jmask == U(i,1));
    gmsens_parcels(idx) = repmat(U(i,5),length(idx),1);
end

npar = length(find(U(:,5)));

% Plot parcels to which we are sensitive
cmap_parc = distinguishable_colors(2500);
cmap_parc(length(cmap_parc)/2+1,:) = 0.9; % to make the rest of the brain grey

fig1 = figure; set(fig1,'units', 'normalized', 'outerposition', [0 0 0.8 0.8], 'Color', [1 1 1]);
subplot(121)
DOTHUB_plotSurfaceImage (gmSurfaceMesh, gmsens_parcels, [0,0], [], cmap_parc)
colorbar off
subplot(122)
DOTHUB_plotSurfaceImage (gmSurfaceMesh, gmsens_parcels, [180,30], [], cmap_parc)
colorbar off
title(['Z axis threshold: ' num2str(z_val) ', Sensitivity threshold: ' num2str(0.1) ', Proportion of nodes: ' num2str(node_prop)])

% Remove Limbic parcel from U and gmsens_parcels (n = 288)
% From bb_parcels_yeo we know that it is only one parcel, so remove it
U(U(:,1)==288, 5) = 0;
gmsens_parcels(gmsens_parcels==288) = 0;
U(U(:,1)==287, 5) = 0;
gmsens_parcels(gmsens_parcels==287) = 0;
npar = length(find(U(:,5)));

fig1 = figure; set(fig1,'units', 'normalized', 'outerposition', [0 0 0.8 0.8], 'Color', [1 1 1]);
subplot(121)
DOTHUB_plotSurfaceImage (gmSurfaceMesh, gmsens_parcels, [0,0], [], cmap_parc)
colorbar off
subplot(122)
DOTHUB_plotSurfaceImage (gmSurfaceMesh, gmsens_parcels, [180,30], [], cmap_parc)
colorbar off
title(['Z axis threshold: ' num2str(z_val) ', Sensitivity threshold: ' num2str(0.1) ', Proportion of nodes: ' num2str(node_prop)])

% Create group level FC matrices for connICA and groupICA
% Select data for analysis
path_dotimg = ('/Volumes/CAM_data/neoLAB/dotimg/ssr');
cd(path_dotimg);
files_dotimg = dir('*.dotimg');

groupFC = zeros(2*npar, 2*npar, length(files_dotimg));
groupFC_z = zeros(2*npar, 2*npar, length(files_dotimg));
groupICA = [];

clc
for nsub = 1:length(files_dotimg)
    nsub
    cd(path_dotimg);
    data = load(files_dotimg(nsub).name, '-mat');
    hbo = data.hbo;
    hbr = data.hbr;
    
    % Computes the mean across nodes to create the parcels
    [hbo_par, hbr_par] = ParcellateNodes(data, gmsens_parcels);

    % Parcel 1 corresponds to all zeros
    hbo_par(:,1) = [];
    hbr_par(:,1) = [];
    
    % FC matrices for connICA
    groupFC(:,:,nsub) = corr([hbo_par, hbr_par]);
    groupFC_z(:,:,nsub) = atanh(groupFC(:,:,nsub));
    
    % Data for groupICA
    % Standardize first
    
    % STANDARDIZATION OXY (mean 0 - variance 1)
    a = repmat(mean(hbo_par,1),[length(hbo_par) 1]);
    b = repmat(std(hbo_par), [length(hbo_par) 1]);
    z_hbo = (hbo_par-a)./b;
    
    % STANDARDIZATION DEOXY (mean 0 - variance 1)
    c = repmat(mean(hbr_par,1),[length(hbr_par) 1]);
    d = repmat(std(hbr_par), [length(hbr_par) 1]);
    z_hbr = (hbr_par-c)./d;
    
    % Save data for each subject
    sub_data.hbo = hbo;
    sub_data.hbr = hbr;
    sub_data.hbo_parc = hbo_par;
    sub_data.hbr_parc = hbr_par;
    sub_data.hbo_zparc = z_hbo;
    sub_data.hbr_zparc = z_hbr;
    sub_data.FC = corr([hbo.gm(:, find(gmsens_parcels)), hbr.gm(:, find(gmsens_parcels))]);    
    sub_data.FC_parc = corr([hbo_par, hbr_par]);
    sub_data.FC_zparc = atanh(corr([hbo_par, hbr_par]));

    save([files_dotimg(nsub).name(1:end-7) '.mat'],'sub_data')

    groupICA = cat(1, [z_hbo, z_hbr], groupICA);
    
end

cd('/Volumes/CAM_data/neoLAB/group')
save('groupICA_allsub_ssr.mat','groupICA')
save('connICA_allsub_ssr.mat','groupFC', 'groupFC_z')


FC_group = tanh(mean(groupFC_z, 3));
imagesc(FC_group, [-1 1]); colormap jet






