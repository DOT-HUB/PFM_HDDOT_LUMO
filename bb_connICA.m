%% =========================
% ======== connICA  ========
% ==========================
% BB 10/20/2020 - b.blanco@bcbl.eu
% Adapted from Amico et al., 2018 - Neuroimage
clear; close all; clc

% Add paths (scripts and toolboxes)
addpath(genpath('/Users/borjablanco/Documents/GitHub/Hyper24RSFC_scripts/scripts/BB_scripts/FastICA_25'));
addpath(genpath('/Users/borjablanco/Documents/GitHub/Hyper24RSFC_scripts/scripts/BB_scripts/icasso122'));

% Load group FC matrices for HbO and HbR
% Dimensions are brain regions x brain regions x subjects
load('/Volumes/CAM_data/neoLAB/group/connICA_allsub_ssr.mat')
FC_hbo = groupFC(1:size(groupFC,1)/2, 1:size(groupFC,1)/2,:);
FC_hbr = groupFC(size(groupFC,1)/2+1:end, size(groupFC,1)/2+1:end,:);

% ConnICA params
configs.numRegions = size(FC_hbo,1); % number of nodes (voxels, channels)
configs.mask = triu(true(configs.numRegions,configs.numRegions),1); % upper triangular mask (symmetric)
configs.numConn = size(FC_hbo,3); % number of connectomes (subjects or sessions)
configs.numEdges = nnz(configs.mask);

oxy_matrix = zeros(configs.numConn,configs.numEdges);
deoxy_matrix = zeros(configs.numConn,configs.numEdges);

% Vectorize each subject's (HbO and HbR) connectomes -> connICA matrix
for i = 1:configs.numConn
    
    aux1 = FC_hbo(:,:,i);
    aux2 = FC_hbr(:,:,i);
    
    oxy_matrix(i,:) = aux1(configs.mask);
    deoxy_matrix(i,:) = aux2(configs.mask);
    
end

% First half columns HbO FC values, second half HbR FC values
% Read hybrid connICA paper - Amico and Go?i, 2018
hybrid_corr = [oxy_matrix deoxy_matrix]';

% Dimensionality reduction (group PCA)
% If the number of samples (ch) is less than or equal to the number of
% components, the n-th principal component will be constant zero (eigenvalue = 0)
[c , s, latent, ~, explained] = pca(hybrid_corr);

% find number of components for different variance thresholds
keep_var = [0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1];

% Compute explained variance across PCA components
% this is the same as the outcome provided by explained in pca
exp_var = cumsum(latent)./sum(latent);

plot(exp_var);figure(gcf);
%variance1 = diff(cumsum(latent)./sum(latent))*100;
%variance1 = cumsum(explained);
M = 100; % can change this number, it's the bumer or rounds of fast ICA to run

% ---------------------- ICASSO PCA selection -------------------
% Compute number of PCA components with the selected explained variance
nComp = zeros(1, length(keep_var));
for i = 1:length(keep_var)
    nComp(i) = length(find(exp_var<keep_var(i)))+1;
end

iq_tot = zeros (length(nComp), max(nComp)); % ICASSO stability idx
hb_cor = zeros (length(nComp), max(nComp)); % correlation between HbO and HbR components
ssd_pca = zeros (length(nComp), max(nComp)); % explained variance pca
ssd_orig = zeros (length(nComp), max(nComp)); % explained variance original
IC_comp = [];
IC_weights = [];

for i = 1:length(nComp)
    clc
    disp(['Number of components retained = ' num2str(nComp(i))]);
    disp(['% of variance retained =  ' num2str(exp_var(nComp(i))*100)]);
    
    % Run PCA
    [coeffs, score, latent] = pca(hybrid_corr,'NumComponents',nComp(i));
    PCA_matrix = score * coeffs'; % PCA reconstructed demeaned data
    PCA_matrix = bsxfun(@plus, PCA_matrix, mean(hybrid_corr, 1)); % plug the mean back
    group_ICA = PCA_matrix;
    
    % Run ICASSO and store results
    % Compute ICASSO stability index (Iq)
    sR = icassoEst('randinit', group_ICA', M,'numOfIC', nComp(i), 'verbose', 'off', 'maxNumIterations', 500, ...
        'approach', 'defl');
    
    % sR contains now estimates. Next:
    %- dissimilarity measure between them is formed
    %- estimates are clustered
    %- a projection for the visualization is computed
    
    % Default similarity measure, clustering and projection:
    % This is the part that takes more time
    sR = icassoExp(sR);
    
    % Visualization & returning results
    icassoShow(sR,'L',nComp(i),'estimate','demixing');
    [iq,ica_weights,w,ica_maps] = icassoShow(sR,'L',nComp(i));    
    %close all   
        
    % Order components by stability index (iq)
    [iq_order, idx_order] = sort(iq, 'descend');
    ica_weights = ica_weights(:, idx_order);
    ica_maps = ica_maps(idx_order, :);
    
    % Store components and iq
    iq_tot(i,1:nComp(i)) = iq_order;
    IC_comp = [IC_comp ica_maps'];
    IC_weights = [IC_weights ica_weights];
    % Calculate and store ssd of each component (Ashby fMRI book)
    % with respect to original data and with respect to data after pca
    for j = 1:nComp(i)
        
        ica_temp = ica_maps;
        ica_temp(j,:) = 0;
        data_rec = ica_weights*ica_temp;
        
        diff_pca = group_ICA - data_rec';
        diff_orig = hybrid_corr - data_rec';
        
        ssd_pca(i,j) = sum(diff_pca(:).^2);
        ssd_orig(i,j) = sum(diff_orig(:).^2);       
    end
    
    ssd_pca(i, :) = ssd_pca(i,:)./sum(ssd_pca(i,:))*100;
    ssd_orig(i,:) = ssd_orig(i,:)./sum(ssd_orig(i,:))*100*keep_var(i);
        
    % Store corr between oxy-deoxy maps
    TC_oxy = ica_maps(:, 1:size(hybrid_corr,1)/2)';
    TC_deoxy = ica_maps(:,(size(hybrid_corr,1)/2)+1:size(hybrid_corr,1))';
    
    temp_cor = corr([TC_oxy, TC_deoxy]);
    temp_diag = zeros(1,nComp(i));
    
    % Select correlation values between equivalent HbO and HbR traits
    for d = 1:nComp(i)   
        temp_diag(1,d) = temp_cor(d, nComp(i)+d);       
    end
    
    hb_cor(i, 1:nComp(i)) = temp_diag;
    
    % Save figure corr between HbO-HbR comp maps
    %figure; imagesc(temp_cor, [-1 1]); colormap jet
    
    %saveas(gcf, ['HbOHbR_connICA_' num2str(keep_var(i)*100) '%.tiff'], 'tiff')
    
    close all
end

% Create and save data structure
% Select correlation values between equivalent HbO and HbR traits
for d = 1:nComp(i)
    temp_diag(1,d) = temp_cor(d, nComp(i)+d);
end

hb_cor(i, 1:nComp(i)) = temp_diag;

% Save figure corr between HbO-HbR comp maps
figure; imagesc(temp_cor, [-1 1]); colormap jet % dia
%saveas(gcf, ['HbOHbR_connICA_' num2str(keep_var(i)*100) '%.tiff'], 'tiff')
close all

% Create and save data structure
connICA.iq_tot = iq_tot;
connICA.hb_cor = hb_cor;
connICA.ssd_pca = ssd_pca;
connICA.ssd_orig = ssd_orig;
connICA.IC_comp = IC_comp;
connICA.weights = IC_weights;
connICA.nComps = nComp;

save('connICA_results_ssr.mat','connICA')


% Create volcano plots
figure;
fig1=gcf;
set(fig1, 'units', 'normalized', 'outerposition', [0 0 1 1], 'Color', [1 1 1]);
subplot(221); imagesc(connICA.iq_tot, [0 1]);box off
set(gca, 'YTick', 1:9 ,'YTickLabels', {'60%','65%','70%','75%','80%','85%','90%','95%', '100%'}, 'FontSize', 20);
ytickangle(30)
subplot(222); imagesc(connICA.hb_cor, [-1 1]); box off
set(gca, 'YTick', 1:9 ,'YTickLabels', {'60%','65%','70%','75%','80%','85%','90%','95%', '100%'}, 'FontSize', 20);
ytickangle(30)
mask_plots = connICA.ssd_orig~=0;
subplot(223); imagesc(connICA.ssd_pca.*mask_plots, [0 40]); box off
set(gca, 'YTick', 1:9 ,'YTickLabels', {'60%','65%','70%','75%','80%','85%','90%','95%', '100%'}, 'FontSize', 20);
ytickangle(30)
subplot(224); imagesc(cumsum(connICA.ssd_orig,2).*mask_plots, [0 100]); box off
set(gca, 'YTick', 1:9 ,'YTickLabels', {'60%','65%','70%','75%','80%','85%','90%','95%', '100%'}, 'FontSize', 20);
ytickangle(30)
colormap jet


figure
imagesc(connICA.weights)
colormap jet








