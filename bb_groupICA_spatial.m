%% ===========================
% ======== GROUP ICA  ========
% ============================
clear; close all; clc

% load('/Volumes/CAM_data/neoLAB/group/parcel_info.mat')
% load('/Users/borjablanco/Documents/GitHub/DOT-HUB_toolbox/Utilities/greyJet.mat')
% origMeshFileName = '/Volumes/CAM_data/neolab/Rob_Mesh/Robmesh_bb.mshs';
% load(origMeshFileName, '-mat')
% % Load parcel info to assign the correct labels to parcels
% load('/Volumes/CAM_data/neoLAB/group/parcel_info.mat')
% load('/Volumes/CAM_data/neoLAB/group/gmSens_mask.mat')

% Add paths (scripts and toolboxes)
addpath(genpath('/Users/borjablanco/Documents/GitHub/Hyper24RSFC_scripts/scripts/BB_scripts/FastICA_25'));
addpath(genpath('/Users/borjablanco/Documents/GitHub/Hyper24RSFC_scripts/scripts/BB_scripts/icasso122'));

% Load concatenated/standardized data
load('/Volumes/CAM_data/neoLAB/group/groupICA_allsub_ssr.mat')


% ----------------------- TEMPORAL HYBRID GROUP ICA METRICS ----------------------
% PRIOR TO ICA
% Dimensionality reduction (group PCA)
% If the number of samples (ch) is less than or equal to the number of
% components, the n-th principal component will be constant zero (eigenvalue = 0)
hybrid_pca = groupICA; % each ch is a row vector

[~ , ~, latent, ~, explained] = pca(hybrid_pca);

% Find number of components for different thresholds
% Check BB RS paper supplementary material
keep_var = [0.6 0.65 0.7 0.75 0.8];

% Compute explained variance across PCA components
% this is the same as the outcome provided by explained in pca
exp_var = cumsum(latent)./sum(latent);

% Plot explained variance
plot(exp_var);

% ---------------------- ICASSO PCA selection -------------------
M = 100; % number of ICASSO randomizations; change this to 100
nComp = zeros(1, length(keep_var));
% Compute number of components for each % of explained variance 
for i = 1:length(keep_var)
    nComp(i) = length(find(exp_var<keep_var(i)))+1; 
end

nica = 12;
% Initialize variables
nch = size(groupICA,2)/2;
iq_tot = zeros (length(nComp), max(nComp)); % ICASSO stability idx
hb_cor = zeros (length(nComp), max(nComp)); % correlation between HbO and HbR components
ssd_pca = zeros (length(nComp), max(nComp)); % explained variance pca
ssd_orig = zeros (length(nComp), max(nComp)); % explained variance original

iq_tot = zeros (length(nComp), nica); % ICASSO stability idx
hb_cor = zeros (length(nComp), nica); % correlation between HbO and HbR components
ssd_pca = zeros (length(nComp), nica); % explained variance pca
ssd_orig = zeros (length(nComp), nica); % explained variance original

IC_timeC = [];
IC_weights_hbo = []; 
IC_weights_hbr = []; 

% For each value in keep_var (nComp)
for i = 1:length(nComp)
    
    % Run PCA
    [coeffs, score, ~ ] = pca(hybrid_pca,'NumComponents',nComp(i));
    PCA_matrix = score * coeffs'; % PCA reconstructed demeaned data
    PCA_matrix = bsxfun(@plus, PCA_matrix,mean(hybrid_pca,1)); % plug the mean back
    
    % Run Icasso and store results
    % This part is computationally heavy
%     sR = icassoEst('randinit', PCA_matrix, M,'numOfIC', nComp(i), 'verbose', 'off', 'maxNumIterations', 500, ...
%         'approach', 'defl');
    sR = icassoEst('randinit', PCA_matrix, M,'numOfIC', nica, 'verbose', 'off', 'maxNumIterations', 500, ...
        'approach', 'defl');
    % sR contains now estimates. Next:
    %- dissimilarity measure between them is formed
    %- estimates are clustered
    %- a projection for visualization is computed
    sR = icassoExp(sR);
    
    % Visualization & returning results
    %icassoShow(sR,'L',nComp(i),'estimate','demixing');
    %[iq,ica_timeCourses,w,ica_weights] = icassoShow(sR,'L',nComp(i));    
    %close all   
    icassoShow(sR,'L',nica,'estimate','demixing');
    [iq,ica_timeCourses,w,ica_weights] = icassoShow(sR,'L',nica);    
    
    % Order components by iq
    [iq_order, idx_order] = sort(iq, 'descend');
    ica_weights = ica_weights';
    ica_timeCourses = ica_timeCourses';
    
    icaSM = ica_weights(:, idx_order); % ICA spatial maps
    icaTC = ica_timeCourses(idx_order, :); % ICA time courses
    
    % Store components and iq
    %iq_tot(i,1:nComp(i)) = iq_order;
    iq_tot(i,1:nica) = iq_order;
    IC_timeC = [IC_timeC icaTC'];
    IC_weights_hbo = [IC_weights_hbo icaSM(1:nch, :)];
    IC_weights_hbr = [IC_weights_hbr icaSM(nch+1:2*nch, :)];

    % Calculate and store ssd of each component - Ashby book
    % with respect to original data and with respect to data after pca
    %for j = 1:nComp(i)
    for j = 1:nica

        % IMPORTANT! The time course are independent, not maps!
        % reconstruct making zero each of the IC time courses
        ica_temp = icaTC;
        ica_temp(j,:) = 0;
        data_rec = icaSM*ica_temp;
        
        % back projection  method 2 - A. Delorme
        %ica_temp = icaSM(:,j)*icaTC(j,:);
        %data_rec = groupICA - ica_temp;
        
        % Compute sum of squared differences
        diff_pca =  PCA_matrix' - data_rec;
        diff_orig = hybrid_pca' - data_rec;
        
        ssd_pca(i,j) = sum(diff_pca(:).^2);
        ssd_orig(i,j) = sum(diff_orig(:).^2);
        
    end  
    
    % Convert it to percentages
    ssd_pca(i, :) = ssd_pca(i,:)./sum(ssd_pca(i,:))*100;
    ssd_orig(i,:) = ssd_orig(i,:)./sum(ssd_orig(i,:))*100*keep_var(i);
        
    % Store negative corr between HbO-HbR spatial maps (expected negative and high)    
    SM_oxy = icaSM (1:nch, :);
    SM_deoxy = icaSM (nch+1:2*nch, :);
    temp_cor = corr(SM_oxy, SM_deoxy);
    
    % we are interested in the main diagonal hbo-hbr but same component
    % correlation should be negative
    %hb_cor(i, 1:nComp(i)) = diag(temp_cor);
    hb_cor(i, 1:nica) = diag(temp_cor);
    % Save figure corr between comp maps
    figure; imagesc(temp_cor, [-1 1]); colormap jet
    
    %saveas(gcf, ['HbOHbR_groupICA_' num2str(keep_var(i)*100) '%.tiff'], 'tiff')
    
    close all
end

% Save data structure
% groupICA_results.iq_tot = iq_tot(1:length(nComp),1:nComp(end));
% groupICA_results.hb_cor = hb_cor(1:length(nComp),1:nComp(end)); % hbcor is the diagonal in the figure above 
% groupICA_results.ssd_pca = ssd_pca(1:length(nComp),1:nComp(end));
% groupICA_results.ssd_orig = ssd_orig(1:length(nComp),1:nComp(end));
% groupICA_results.IC_comp = IC_timeC;
% groupICA_results.IC_weights_hbo = IC_weights_hbo;
% groupICA_results.IC_weights_hbr = IC_weights_hbr;
% groupICA_results.nComp = nComp;

groupICA_results.iq_tot = iq_tot(1:length(nComp),1:nica);
groupICA_results.hb_cor = hb_cor(1:length(nComp),1:nica); % hbcor is the diagonal in the figure above 
groupICA_results.ssd_pca = ssd_pca(1:length(nComp),1:nica);
groupICA_results.ssd_orig = ssd_orig(1:length(nComp),1:nica);
groupICA_results.IC_comp = IC_timeC;
groupICA_results.IC_weights_hbo = IC_weights_hbo;
groupICA_results.IC_weights_hbr = IC_weights_hbr;
groupICA_results.nica = nica;

save('groupICA_spatial_results_ssr_nica.mat','groupICA_results')

% Create volcano plots
figure;
fig1=gcf;
set(fig1, 'units', 'normalized', 'outerposition', [0 0 1 1], 'Color', [1 1 1]);
subplot(221); imagesc(groupICA_results.iq_tot, [0 1]); box off
set(gca, 'YTick', 1:8 ,'YTickLabels', {'60%','65%','70%','75%','80%','85%','90%','95%'}, 'FontSize', 20);
ytickangle(30)
subplot(222); imagesc(groupICA_results.hb_cor, [-1 0]); box off
set(gca, 'YTick', 1:8 ,'YTickLabels', {'60%','65%','70%','75%','80%','85%','90%','95%'}, 'FontSize', 20);
ytickangle(30)
mask_plots = groupICA_results.ssd_orig~=0;
subplot(223); imagesc(groupICA_results.ssd_pca.*mask_plots, [0 10]); box off
set(gca, 'YTick', 1:8 ,'YTickLabels', {'60%','65%','70%','75%','80%','85%','90%','95%'}, 'FontSize', 20);
ytickangle(30)
subplot(224); imagesc(cumsum(groupICA_results.ssd_orig,2).*mask_plots, [0 100]); box off
set(gca, 'YTick', 1:8 ,'YTickLabels', {'60%','65%','70%','75%','80%','85%','90%','95%'}, 'FontSize', 20);
ytickangle(30)
colormap jet

figure; imagesc(corr(groupICA_results.IC_comp), [-1 1]);
colormap jet
