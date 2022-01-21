% Dual regression and internetwork connectivity
clear; close all; clc
path_data ='/Volumes/CAM_data/neoLAB/group/Data4DualRegression';
cd(path_data)
load('/Volumes/CAM_data/neoLAB/group/groupICA_spatial_results_ssr_nica.mat')

% Extract data
% Convert weights of icaSM to z-scores for plots
sub_hbo = dir('*hbo*');
nsub = length(sub_hbo);
%cumsum(groupICA_results.nComp) -> 107 = the first 5 PCA trhesholds (80%)
icaSM_hbo = groupICA_results.IC_weights_hbo(:, 1:60);
icaSM_std_hbo = zscore(icaSM_hbo);

sub_hbr = dir('*hbr*');
%cumsum(groupICA_results.nComp) -> 107 = the first 5 PCA trhesholds (80%)
icaSM_hbr = groupICA_results.IC_weights_hbr(:, 1:60);
icaSM_std_hbr = zscore(icaSM_hbr);

% Initialize variables
ch = size(icaSM_hbo,1); % n parcels
max_dur = 947; % to make all subjects the same duration
nComp = size(icaSM_std_hbo, 2);

group_sm_hbo = zeros(nComp, ch, nsub);
group_tc_hbo = zeros(max_dur, nComp, nsub);
group_sm_hbr = zeros(nComp, ch, nsub);
group_tc_hbr = zeros(max_dur, nComp, nsub);
% Remove components with low Iq or correlation
%icaSM_std(:,2) = []; 
clc

for i = 1:nsub
    
    load(sub_hbo(i).name)
    load(sub_hbr(i).name)

    % The model is what determines the orientation of the matrices
    % because it has only one direction in which it is correct.
    % METHOD 1 - Two step regression
    % DATA input Y = preprocessed data from one subject (ch x time points)
    % MODEL input X = ICA Spatial Maps from group ICA (ch x ICs)
    % BETAS output = IC associated time courses for each subject (ICs x time points)   
    Y_hbo = hbo(1:max_dur,:)'; % all same length for storing
    X_hbo = icaSM_std_hbo;
    sub_tc_hbo = pinv(X_hbo)*Y_hbo;
    
    Y_hbr = hbr(1:max_dur,:)'; % all same length for storing
    X_hbr = icaSM_std_hbr;
    sub_tc_hbr = pinv(X_hbr)*Y_hbr;
    
    % DATA input = preprocessed data from one subject (time points x ch)
    % MODEL input = IC associated time courses for each subject (time points x ICs)
    % BETAS output = IC associated spatial maps for each subject (ICs x ch)
    Y_hbo = Y_hbo';
    X_hbo = sub_tc_hbo';
    sub_sm_hbo = pinv(X_hbo)*Y_hbo;
    
    Y_hbr = Y_hbr';
    X_hbr = sub_tc_hbr';
    sub_sm_hbr = pinv(X_hbr)*Y_hbr;

    % Store results
    group_sm_hbo(:,:,i) = sub_sm_hbo;
    group_tc_hbo(:,:,i) = sub_tc_hbo';
    
    group_sm_hbr(:,:,i) = sub_sm_hbr;
    group_tc_hbr(:,:,i) = sub_tc_hbr';
    
        
    % METHOD 2 - The same but based on (Erhardt, 2011) paper
    % sub_tc =  sub_std3(:,:,i) * icaSM_std' * inv((icaSM_std*icaSM_std'));
    % sub_sm = inv(sub_tc'*sub_tc)*sub_tc' * sub_std3(:,:,i);
    
    % METHOD 3 - Same outcome
    % sub_sm = linsolve(sub_tc, sub_std3(:,:,i));
    
end

% INTERNETWORK CONNECTIVITY
figure
mat_all = zeros(21,21,14);
for i = 1:size(group_sm_hbo,3)
    subplot(4,4,i)
    
    imagesc(corr(group_tc_hbo(:,32:52,i)), [-1 1])
    
    mat_all(:,:,i) = corr(group_tc_hbo(:,32:52,i));
    colormap jet

end

figure
imagesc(tanh(mean(atanh(mat_all),3)), [-1 1])
colormap jet
set(gca,'YTick',1:21, 'YtickLabel', {'occ ext RH', 'occ ext2 LH', 'frontal RH'...
    'occ ext RH', 'occ ext2 RH', 'frontal LH', 'visC 2', 'visP', 'control RH'...
    'x', 'x', 'occ sup RH', 'salience + DMN', 'x', 'salience LH', 'control LH'...
    'visC', 'occ sup LH', 'x', 'salience RH', 'control + DMN'},'XTick',1:21, 'XTickLabel', ...
    {'occ ext RH', 'occ ext2 LH', 'frontal RH'...
    'occ ext RH', 'occ ext2 RH', 'frontal LH', 'visC 2', 'visP', 'control RH'...
    'x', 'x', 'occ sup RH', 'salience + DMN', 'x', 'salience LH', 'control LH'...
    'visC', 'occ sup LH', 'x', 'salience RH', 'control + DMN'}, 'XTickLabelRotation', 25)







