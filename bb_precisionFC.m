
% Script for precision functional mapping on HD-DOT data
% 
%
% Created by Borja Blanco 07/2021
% bb579@cam.ac.uk

% Set up environment
clear; close all; clc
addpath(genpath('/Users/borjablanco/Documents/GitHub/Hyper24RSFC_scripts/scripts/BB_scripts'))
path_data = '/Volumes/CAM_data/neoLAB/dotimg/ssr';
cd(path_data)
sub = dir('*ssr.mat');

% Load parcel info
load('/Volumes/CAM_data/neoLAB/group/gmSens_parcels.mat');

% Load Schaefer labels
fid = fopen('/Volumes/CAM_data/neoLAB/group/Schaefer2018_1000Parcels_17Networks_order.txt');
label_info = textscan(fid,'%f%s%f%f%f%f');
fclose(fid);

% Assign labels to parcels
gm_parc = unique(gmsens_parcels); % find unique parcels
gm_parc(gm_parc==0) = []; % Remove parcel of 0s

parc_idx = [];
for i = 1:length(gm_parc)
   parc_idx(i) = find(label_info{:,1}==gm_parc(i));
end

% Assign Yeo labels to parcels
parc_label = label_info{1,2};
parc_label = parc_label(parc_idx);

% Rename labels to make them shorter and unique
for i = 1:length(parc_label)
    temp_label = parc_label{i,1};
    temp_label = temp_label(1:18);
    parc_label{i,1} = temp_label;
end

% Parcels correspond to 12 Yeo networks (6 LH and 6 RH) 
[a,b, net_labels] = unique(parc_label);
parc_idx = parc_idx';
parc_idx(:,2) = net_labels;

% Merge LH and RH networks for analysis
[a,b,c]=unique(parc_idx(:,2), 'stable');
a_idx = a;
a_idx(length(a)/2+1:end) = a(1:length(a)/2);
parc_net = parc_idx;

for i = 1:length(a)
   parc_net(parc_net(:,2) == a_idx(i)+length(a)/2,2) = a_idx(i);
end

% ---------------------------------------------------------
%  --------------------- PRECISION FC
% PART 1 - Compute correlation between FC matrices of each session
% and correlation with ground truth matrix, which is
% the FC matrix based on concatenated data from all sessions
% This is done for each network independently in the analysis
hbo_group = [];
hbr_group = [];
cd(path_data)
max_size = 947; % some participants have 947 and others 948 samples
for nsub = 1:length(sub)
    nsub
    load(sub(nsub).name, '-mat');
    
    % Load parcel data normalized to mean 0 and std 1
    hbo_sub = sub_data.hbo_zparc(1:max_size,:);
    hbr_sub = sub_data.hbr_zparc(1:max_size,:);
    
    % Store data for each participant in group level matrix
    hbo_group = cat(3,hbo_group, hbo_sub);
    hbr_group = cat(3,hbr_group, hbr_sub);
    clc
end

% Already computed, takes time, load variables for speed
load('/Volumes/CAM_data/neoLAB/group/precisionFC_group_data_ssr.mat')

% Initialize variables to store data at the individual and group level
visC_hbo = []; visC_hbr = [];
visP_hbo = []; visP_hbr = [];
dors_hbo = []; dors_hbr = [];
sal_hbo = []; sal_hbr = [];
cont_hbo = []; cont_hbr = [];
defa_hbo = []; defa_hbr = [];
all_hbo  = []; all_hbr = [];

group_visC_hbo = []; group_visC_hbr = [];
group_visP_hbo = []; group_visP_hbr = [];
group_dors_hbo = []; group_dors_hbr = [];
group_sal_hbo = []; group_sal_hbr = [];
group_cont_hbo = []; group_cont_hbr = [];
group_defa_hbo = []; group_defa_hbr = [];
group_all_hbo = []; group_all_hbr = [];

% Separate parcels into networks based on Schaefer - Yeo
% and create group matrices for networkwise ground truth
for nsub = 1:length(sub)
    nsub
    visC_hbo (:, :, nsub) = hbo_group(:, parc_net(:,2)==a_idx(1), nsub);
    visC_hbr (:, :, nsub) = hbr_group(:, parc_net(:,2)==a_idx(1), nsub);
    
    visP_hbo (:, :, nsub) = hbo_group(:, parc_net(:,2)==a_idx(2), nsub);
    visP_hbr (:, :, nsub) = hbr_group(:, parc_net(:,2)==a_idx(2), nsub);
    
    dors_hbo (:, :, nsub) = hbo_group(:, parc_net(:,2)==a_idx(3), nsub);
    dors_hbr (:, :, nsub) = hbr_group(:, parc_net(:,2)==a_idx(3), nsub);
    
    sal_hbo (:, :, nsub) = hbo_group(:, parc_net(:,2)==a_idx(4), nsub);
    sal_hbr (:, :, nsub) = hbr_group(:, parc_net(:,2)==a_idx(4), nsub);
    
    cont_hbo (:, :, nsub) = hbo_group(:, parc_net(:,2)==a_idx(5), nsub);
    cont_hbr (:, :, nsub) = hbr_group(:, parc_net(:,2)==a_idx(5), nsub);
    
    defa_hbo (:, :, nsub) = hbo_group(:, parc_net(:,2)==a_idx(6), nsub);
    defa_hbr (:, :, nsub) = hbr_group(:, parc_net(:,2)==a_idx(6), nsub);
          
    all_hbo (:, :, nsub) = hbo_group(:, :, nsub);
    all_hbr (:, :, nsub) = hbr_group(:, :, nsub);
    
    group_visC_hbo = [group_visC_hbo; squeeze(visC_hbo(:,:,nsub))];
    group_visC_hbr = [group_visC_hbr; squeeze(visC_hbr(:,:,nsub))];
    
    group_visP_hbo = [group_visP_hbo; squeeze(visP_hbo(:,:,nsub))];
    group_visP_hbr = [group_visP_hbr; squeeze(visP_hbr(:,:,nsub))];
    
    group_dors_hbo = [group_dors_hbo; squeeze(dors_hbo(:,:,nsub))];
    group_dors_hbr = [group_dors_hbr; squeeze(dors_hbr(:,:,nsub))];
    
    group_sal_hbo = [group_sal_hbo; squeeze(sal_hbo(:,:,nsub))];
    group_sal_hbr = [group_sal_hbr; squeeze(sal_hbr(:,:,nsub))];
    
    group_cont_hbo = [group_cont_hbo; squeeze(cont_hbo(:,:,nsub))];
    group_cont_hbr = [group_cont_hbr; squeeze(cont_hbr(:,:,nsub))];
    
    group_defa_hbo = [group_defa_hbo; squeeze(defa_hbo(:,:,nsub))];
    group_defa_hbr = [group_defa_hbr; squeeze(defa_hbr(:,:,nsub))];
    
    group_all_hbo = [group_all_hbo; squeeze(all_hbo(:,:,nsub))];
    group_all_hbr = [group_all_hbr; squeeze(all_hbr(:,:,nsub))];
end

% Compute and plot SIMILARITY MATRICES
% for all sessions and all parcels in each network
visC_hbo_plots = reshape(visC_hbo, [size(hbo_group, 1), size(visC_hbo,2)*length(sub)]);
visP_hbo_plots = reshape(visP_hbo, [size(hbo_group, 1), size(visP_hbo,2)*length(sub)]);
dors_hbo_plots = reshape(dors_hbo, [size(hbo_group, 1), size(dors_hbo,2)*length(sub)]);
sal_hbo_plots = reshape(sal_hbo, [size(hbo_group, 1), size(sal_hbo,2)*length(sub)]);
cont_hbo_plots = reshape(cont_hbo, [size(hbo_group, 1), size(cont_hbo,2)*length(sub)]);
defa_hbo_plots = reshape(defa_hbo, [size(hbo_group, 1), size(defa_hbo,2)*length(sub)]);
all_hbo_plots = reshape(all_hbo, [size(hbo_group, 1), size(all_hbo,2)*length(sub)]);

visC_hbr_plots = reshape(visC_hbr, [size(hbr_group, 1), size(visC_hbr,2)*length(sub)]);
visP_hbr_plots = reshape(visP_hbr, [size(hbr_group, 1), size(visP_hbr,2)*length(sub)]);
dors_hbr_plots = reshape(dors_hbr, [size(hbr_group, 1), size(dors_hbr,2)*length(sub)]);
sal_hbr_plots = reshape(sal_hbr, [size(hbr_group, 1), size(sal_hbr,2)*length(sub)]);
cont_hbr_plots = reshape(cont_hbr, [size(hbr_group, 1), size(cont_hbr,2)*length(sub)]);
defa_hbr_plots = reshape(defa_hbr, [size(hbr_group, 1), size(defa_hbr,2)*length(sub)]);
all_hbr_plots = reshape(all_hbr, [size(hbr_group, 1), size(all_hbr,2)*length(sub)]);

% Plot average within and between session similarity
bb_similarity_plot(visC_hbo_plots, visC_hbr_plots, size(visC_hbo,2), length(sub), 'Visual Central')
bb_similarity_plot(visP_hbo_plots, visP_hbr_plots, size(visP_hbo,2), length(sub), 'Visual Peripheral')
bb_similarity_plot(dors_hbo_plots, dors_hbr_plots, size(dors_hbo,2), length(sub), 'Dorsal Attention')
bb_similarity_plot(sal_hbo_plots, sal_hbr_plots, size(sal_hbo,2), length(sub), 'Salience')
bb_similarity_plot(cont_hbo_plots, cont_hbr_plots, size(cont_hbo,2), length(sub), 'Control')
bb_similarity_plot(defa_hbo_plots, defa_hbr_plots, size(defa_hbo,2), length(sub), 'Default-Mode Network')
bb_similarity_plot(all_hbo_plots, all_hbr_plots, size(all_hbo,2), length(sub), 'All')

% possible colormaps: lapaz, nuuk, tokyo, bilbao
% Similarity across sessions, and with respect to ground truth
[vc_hbo, vc_hbr] = bb_similarity(visC_hbo, group_visC_hbo, visC_hbr, group_visC_hbr, 'Visual central');
[vp_hbo, vp_hbr] = bb_similarity(visP_hbo, group_visP_hbo, visP_hbr, group_visP_hbr, 'Visual peripheral');
[do_hbo, do_hbr] = bb_similarity(dors_hbo, group_dors_hbo, dors_hbr, group_dors_hbr, 'Dorsal attention');
[sa_hbo, sa_hbr] = bb_similarity(sal_hbo, group_sal_hbo, sal_hbr, group_sal_hbr, 'Salience');
[co_hbo, co_hbr] = bb_similarity(cont_hbo, group_cont_hbo, cont_hbr, group_cont_hbr, 'Control');
[de_hbo, de_hbr] = bb_similarity(defa_hbo, group_defa_hbo, defa_hbr, group_defa_hbr, 'Default-mode');
[ap_hbo, ap_hbr] = bb_similarity(all_hbo, group_all_hbo, all_hbr, group_all_hbr, 'All parcels');

% SECOND PART
% STABILITY FOR DIFFERENT RECORDING DURATIONS
% Compute ground truth matrix (all sessions, except the one assessed)
bb_stability(visC_hbo, group_visC_hbo, visC_hbr, group_visC_hbr, 'Visual central')
bb_stability(visP_hbo, group_visP_hbo, visP_hbr, group_visP_hbr, 'Visual peripheral')
bb_stability(dors_hbo, group_dors_hbo, dors_hbr, group_dors_hbr, 'Dorsal attention')
bb_stability(sal_hbo, group_sal_hbo, sal_hbr, group_sal_hbr, 'Salience')
bb_stability(cont_hbo, group_cont_hbo, cont_hbr, group_cont_hbr, 'Control')
bb_stability(defa_hbo, group_defa_hbo, defa_hbr, group_defa_hbr, 'Default-mode')
bb_stability(all_hbo, group_all_hbo, all_hbr, group_all_hbr, 'All parcels')

% % Violin plots across sessions
% s_hbo = [];
% s_hbr = [];
% 
% % Extract similarity of each session with the ground truth
% for i = 1:14
%     s_hbo(:,i) = [vc_hbo(i,15), vp_hbo(i,15), do_hbo(i,15), sa_hbo(i,15), co_hbo(i,15), de_hbo(i,15)];
%     s_hbr(:,i) = [vc_hbr(i,15), vp_hbr(i,15), do_hbr(i,15), sa_hbr(i,15), co_hbr(i,15), de_hbr(i,15)];
% end
% 
% figure
% subplot(121)
% boxplot(s_hbo')
% title('HbO')
% xlabel('Networks')
% ylabel('Similarity (r)')
% subplot(122)
% boxplot(s_hbr')
% title('HbR')
% xlabel('Networks')
% ylabel('Similarity (r)')
% 
% % Compute intersession mean values for each network
% vc_inter_hbo = triu(vc_hbo,1); vc_inter_hbo(:,15) = 0; vc_inter_hbo = vc_inter_hbo(vc_inter_hbo~=0);
% vc_inter_hbr = triu(vc_hbr,1); vc_inter_hbr(:,15) = 0; vc_inter_hbr = vc_inter_hbr(vc_inter_hbr~=0);
% vp_inter_hbo = triu(vp_hbo,1); vp_inter_hbo(:,15) = 0; vp_inter_hbo = vp_inter_hbo(vp_inter_hbo~=0);
% vp_inter_hbr = triu(vp_hbr,1); vp_inter_hbr(:,15) = 0; vp_inter_hbr = vp_inter_hbr(vp_inter_hbr~=0);
% do_inter_hbo = triu(do_hbo,1); do_inter_hbo(:,15) = 0; do_inter_hbo = do_inter_hbo(do_inter_hbo~=0);
% do_inter_hbr = triu(do_hbr,1); do_inter_hbr(:,15) = 0; do_inter_hbr = do_inter_hbr(do_inter_hbr~=0);
% sa_inter_hbo = triu(sa_hbo,1); sa_inter_hbo(:,15) = 0; sa_inter_hbo = sa_inter_hbo(sa_inter_hbo~=0);
% sa_inter_hbr = triu(sa_hbr,1); sa_inter_hbr(:,15) = 0; sa_inter_hbr = sa_inter_hbr(sa_inter_hbr~=0);
% co_inter_hbo = triu(co_hbo,1); co_inter_hbo(:,15) = 0; co_inter_hbo = co_inter_hbo(co_inter_hbo~=0);
% co_inter_hbr = triu(co_hbr,1); co_inter_hbr(:,15) = 0; co_inter_hbr = co_inter_hbr(co_inter_hbr~=0);
% de_inter_hbo = triu(de_hbo,1); de_inter_hbo(:,15) = 0; de_inter_hbo = de_inter_hbo(de_inter_hbo~=0);
% de_inter_hbr = triu(de_hbr,1); de_inter_hbr(:,15) = 0; de_inter_hbr = de_inter_hbr(de_inter_hbr~=0);

% GRAPH METRICS
% For debugging ---------------------
load('/Volumes/CAM_data/neoLAB/group/precisionFC_group_data_ssr.mat')
% -----------------------------------

% Define parameters for sliding window
sf = 1.315; % downsampled from 5.26Hz/4 = 1.315
data_size = 936; % samples, to have integer outputs below
% define time windows (around 60 seconds time windows)
% 1,2,3,4,5,6,7,8,9,10,11,12 minutes
tw = 1:12;
twp = data_size/length(tw);
tw_vect = tw.*twp;

% Initialize variables
% 2 graph metrics
% 144 time windows in total
% 50 sparsity thresholds
% 14 sessions
graphmet_hbo = zeros(5,144,50,14);
graphmet_hbr = zeros(5,144,50,14);

parpool('local', 4); % use 4 CPUs
parallelpool = gcp;

fprintf('Computing graph metrics, progress:\n');
parfor_progress(size(hbo_group,3));
parfor nsub = 1:size(hbo_group,3)
    
    % Compute graph metrics for each session (nsub)
    % Based on FC with all the parcels - not by network
    % Reduce data to 936 samples - easier to create time windows
    ses_fc_hbo = hbo_group(1:twp*length(tw),:,nsub)';
    ses_fc_hbr = hbr_group(1:twp*length(tw),:,nsub)';
    
    % Graph metrics computed on 1 to 12 minutes time windows with 30 seconds overlap
    [graphmet_hbo(:,:,:,nsub), graphmet_hbr(:,:,:,nsub), ~] = bb_graphTimeWindow(ses_fc_hbo, ses_fc_hbr, tw_vect);
    
    parfor_progress;
end

parfor_progress(0);

% These variables have been computed
load('/Volumes/CAM_data/neoLAB/group/graph_met1_50_prop.mat')
load('/Volumes/CAM_data/neoLAB/group/graph_met1_50_GroundT.mat')

ge_hbo = squeeze(graphmet_hbo(1,:,:,:));
ge_hbr = squeeze(graphmet_hbr(1,:,:,:));

cl_hbo = squeeze(graphmet_hbo(2,:,:,:));
cl_hbr = squeeze(graphmet_hbr(2,:,:,:));

ml_hbo = squeeze(graphmet_hbo(3,:,:,:));
ml_hbr = squeeze(graphmet_hbr(3,:,:,:));

cc_hbo = squeeze(graphmet_hbo(4,:,:,:));
cc_hbr = squeeze(graphmet_hbr(4,:,:,:));

as_hbo = squeeze(graphmet_hbo(5,:,:,:));
as_hbr = squeeze(graphmet_hbr(5,:,:,:));

% Compute median values in each individual for each time window
med_ge_hbo = [];
med_ge_hbr = [];
med_cl_hbo = [];
med_cl_hbr = [];
med_ml_hbo = [];
med_ml_hbr = [];
med_cc_hbo = [];
med_cc_hbr = [];
med_as_hbo = [];
med_as_hbr = [];

tw = 1:12;
ntw = [];

for n = 1:size(ge_hbo,3)
    for i = 1:length(tw)
        tw_idx = find(tw_order==i);
        med_ge_hbo(:,i,n) = median(ge_hbo(tw_idx,:,n),1);
        med_ge_hbr(:,i,n) = median(ge_hbr(tw_idx,:,n),1);
        med_cl_hbo(:,i,n) = median(cl_hbo(tw_idx,:,n),1);
        med_cl_hbr(:,i,n) = median(cl_hbr(tw_idx,:,n),1);
        med_ml_hbo(:,i,n) = median(ml_hbo(tw_idx,:,n),1);
        med_ml_hbr(:,i,n) = median(ml_hbr(tw_idx,:,n),1);
        med_cc_hbo(:,i,n) = median(cc_hbo(tw_idx,:,n),1);
        med_cc_hbr(:,i,n) = median(cc_hbr(tw_idx,:,n),1);
        med_as_hbo(:,i,n) = median(as_hbo(tw_idx,:,n),1);
        med_as_hbr(:,i,n) = median(as_hbr(tw_idx,:,n),1);
        ntw(i) = length(tw_idx);
    end
end


% ALREADY COMPUTED, DATA LOADED ABOVE

% % Compute graph metrics for ground T data for comparison
% % Compute graph metrics at different density thresholds
% % Create subject specific ground T matrices
% tp = 947;
% ini_ses = 1:tp:size(group_all_hbo,1); 
% end_ses = tp:tp:size(group_all_hbo,1);
% 
% % Create upper triangular mask
% mask = triu(true(size(group_all_hbo,2),size(group_all_hbo,2)),1); 
% 
% groundT_ses_hbo = [];
% groundT_fc_hbo = [];
% groundT_ses_hbr = [];
% groundT_fc_hbr = [];
% 
% for i = 1:size(ge_hbo,3)
%     
%     % Concatenate all sessions except the actual session - JU
%     temp_hbo = group_all_hbo;
%     temp_hbo(ini_ses(i):end_ses(i),:) = [];
%     groundT_ses_hbo(:,:,i) = temp_hbo;
%     % Compute FC matrix of ground T data for each session
%     groundT_fc_hbo (:,:,i) = corr(groundT_ses_hbo(:,:,i));
%     
%     temp_hbr = group_all_hbr;
%     temp_hbr(ini_ses(i):end_ses(i),:) = [];
%     groundT_ses_hbr(:,:,i) = temp_hbr;
%     % Compute FC matrix of ground T data for each session
%     groundT_fc_hbr(:,:,i) = corr(groundT_ses_hbr(:,:,i));
%     
% end
% 
% % Density x Time window x Subjects
% dens = 1:1:50;
% 
% % Initialize variables
% ge_hbo_all = zeros(length(dens), size(groundT_fc_hbo,3));
% ge_hbr_all = zeros(length(dens), size(groundT_fc_hbo,3));
% cl_hbo_all = zeros(length(dens), size(groundT_fc_hbo,3));
% cl_hbr_all = zeros(length(dens), size(groundT_fc_hbo,3));
% ml_hbo_all = zeros(length(dens), size(groundT_fc_hbo,3));
% ml_hbr_all = zeros(length(dens), size(groundT_fc_hbo,3));
% cc_hbo_all = zeros(length(dens), size(groundT_fc_hbo,3));
% cc_hbr_all = zeros(length(dens), size(groundT_fc_hbo,3));
% as_hbo_all = zeros(length(dens), size(groundT_fc_hbo,3));
% as_hbr_all = zeros(length(dens), size(groundT_fc_hbo,3));
% 
% 
% for n = 1:size(groundT_fc_hbo,3)
%     
%     temp_hbo = squeeze(groundT_fc_hbo(:,:,n));
%     temp_hbr = squeeze(groundT_fc_hbr(:,:,n));
%     
%     for th = 1:length(dens)
%         vect_hbo = triu(temp_hbo,1);
%         vect_hbo = sort(vect_hbo(vect_hbo~=0), 'descend');
%         th_hbo = vect_hbo(round(dens(th)/100*length(vect_hbo))+1);
%         
%         vect_hbr = triu(temp_hbr,1);
%         vect_hbr = sort(vect_hbr(vect_hbr~=0), 'descend');
%         th_hbr = vect_hbr(round(dens(th)/100*length(vect_hbr))+1);
%         
%         % Apply threshold to fc matrix
%         ses_hbo = (temp_hbo>=th_hbo).*temp_hbo;
%         ses_hbr = (temp_hbr>=th_hbr).*temp_hbr;
%         
%         % Remove NaN from main diag
%         ses_hbo(isnan(ses_hbo)) = 1;
%         ses_hbr(isnan(ses_hbr)) = 1;
%         
%         % Compute graph metrics
%         % Global efficiency
%         ge_hbo_all (th,n) = efficiency_wei((ses_hbo>0).*ses_hbo,0);
%         ge_hbr_all (th,n) = efficiency_wei((ses_hbr>0).*ses_hbr,0);
%         
%         % Modularity
%         [~, cl_hbo_all(th,n)] = community_louvain(ses_hbo,1, [], 'negative_asym');
%         [~, cl_hbr_all(th,n)] = community_louvain(ses_hbr,1, [], 'negative_asym');
%         
%         [~, ml_hbo_all(th,n)] = community_louvain((ses_hbo>0).*ses_hbo,1, [], 'modularity');
%         [~, ml_hbr_all(th,n)] = community_louvain((ses_hbr>0).*ses_hbr,1, [], 'modularity');
%         
%         % Clustering coefficient
%         cc_hbo_all(th,n) = mean(clustering_coef_wu((ses_hbo>0).*ses_hbo));
%         cc_hbr_all(th,n) = mean(clustering_coef_wu((ses_hbr>0).*ses_hbr));
%         
%         % Assortativity
%         as_hbo_all(th,n) = assortativity_wei((ses_hbo>0).*ses_hbo, 0);
%         as_hbr_all(th,n) = assortativity_wei((ses_hbr>0).*ses_hbr, 0);
%               
%     end
%     clc
%     n
% end
% 
% % Plot figures

% 1 - % Changes with density (understand why to include this - maybe supp.)
for i = 1:50
    figure
    subplot(251)
    plot(squeeze(med_ge_hbo(i,:,:)), 'LineWidth', 1)
    subplot(252)
    plot(squeeze(med_cl_hbo(i,:,:)), 'LineWidth', 1)
    subplot(253)
    plot(squeeze(med_ml_hbo(i,:,:)), 'LineWidth', 1)
    subplot(254)
    plot(squeeze(med_cc_hbo(i,:,:)), 'LineWidth', 1)
    subplot(255)
    plot(squeeze(med_as_hbo(i,:,:)), 'LineWidth', 1)
    
    subplot(256)
    plot(squeeze(med_ge_hbr(i,:,:)), 'LineWidth', 1)
    subplot(257)
    plot(squeeze(med_cl_hbr(i,:,:)), 'LineWidth', 1)
    subplot(258)
    plot(squeeze(med_ml_hbr(i,:,:)), 'LineWidth', 1)
    subplot(259)
    plot(squeeze(med_cc_hbr(i,:,:)), 'LineWidth', 1)
    subplot(2,5,10)
    plot(squeeze(med_as_hbr(i,:,:)), 'LineWidth', 1)
    pause
    close all
end

% 1 - % Changes with time
time = 12;
fig1 = figure; set(fig1, 'units', 'normalized', 'outerposition', [0 0 1 1], 'Color', [1 1 1]);
subplot(231)
plot(squeeze(med_ge_hbo(:,time,:)), 'LineWidth', 1)
ylim([0 1]);xlim([0.5 50.5]); 
xlabel('Graph density (%)')
set(gca, 'Fontsize', 18); box off
ylabel('Global efficiency', 'Fontsize', 20, 'FontWeight', 'Bold')

%subplot(252)
%plot(squeeze(med_cl_hbo(:,time,:)), 'LineWidth', 1)
%ylim([0 1]);xlim([0.5 50.5]); axis square; title('Modularity1')
%xlabel('Density')
subplot(232)
plot(squeeze(med_ml_hbo(:,time,:)), 'LineWidth', 1)
ylim([0 1]);xlim([0.5 50.5]); 
xlabel('Graph density (%)')
set(gca, 'Fontsize', 18); box off
ylabel('Modularity', 'Fontsize', 20, 'FontWeight', 'Bold')

subplot(233)
plot(squeeze(med_cc_hbo(:,time,:)), 'LineWidth', 1)
ylim([0 1]);xlim([0.5 50.5]); 
xlabel('Graph density (%)')
set(gca, 'Fontsize', 18); box off
ylabel('Clustering coefficient', 'Fontsize', 20, 'FontWeight', 'Bold')

%subplot(255)
%plot(squeeze(med_as_hbo(:,time,:)), 'LineWidth', 1)
%ylim([0 1]);xlim([0.5 50.5]); axis square; title('Assortativity')
%xlabel('Density')

subplot(234)
plot(squeeze(med_ge_hbr(:,time,:)), 'LineWidth', 1)
ylim([0 1]);xlim([0.5 50.5]); 
xlabel('Graph density (%)')
set(gca, 'Fontsize', 18); box off
ylabel('Global efficiency', 'Fontsize', 20, 'FontWeight', 'Bold')

subplot(235)
plot(squeeze(med_cl_hbr(:,time,:)), 'LineWidth', 1)
ylim([0 1]);xlim([0.5 50.5]);
xlabel('Graph density (%)')
set(gca, 'Fontsize', 18); box off
ylabel('Modularity', 'Fontsize', 20, 'FontWeight', 'Bold')

%subplot(258)
%plot(squeeze(med_ml_hbr(:,time,:)), 'LineWidth', 1)
%ylim([0 1]);xlim([0.5 50.5]); axis square;
%xlabel('Density')
subplot(236)
plot(squeeze(med_cc_hbr(:,time,:)), 'LineWidth', 1)
ylim([0 1]);xlim([0.5 50.5]); 
xlabel('Graph density (%)')
set(gca, 'Fontsize', 18); box off
ylabel('Clustering coefficient', 'Fontsize', 20, 'FontWeight', 'Bold')

%subplot(2,5,10)
%plot(squeeze(med_as_hbr(:,time,:)), 'LineWidth', 1)
%ylim([0 1]);xlim([0.5 50.5]); axis square;
%xlabel('Density')

% 1b
for i = 1:50-1
    figure
    subplot(251)
    imagesc(corr(squeeze(med_ge_hbo(i,:,:)),squeeze(med_ge_hbo(i+1,:,:))), [-1 1])
    subplot(252)
    imagesc(corr(squeeze(med_cl_hbo(i,:,:)),squeeze(med_cl_hbo(i+1,:,:))), [-1 1])
    subplot(253)
    imagesc(corr(squeeze(med_ml_hbo(i,:,:)),squeeze(med_ml_hbo(i+1,:,:))), [-1 1])
    subplot(254)
    imagesc(corr(squeeze(med_cc_hbo(i,:,:)),squeeze(med_cc_hbo(i+1,:,:))), [-1 1])
    subplot(255)
    imagesc(corr(squeeze(med_as_hbo(i,:,:)),squeeze(med_as_hbo(i+1,:,:))), [-1 1])
    
    subplot(256)
    imagesc(corr(squeeze(med_ge_hbr(i,:,:)),squeeze(med_ge_hbr(i+1,:,:))), [-1 1])
    subplot(257)
    imagesc(corr(squeeze(med_cl_hbr(i,:,:)),squeeze(med_cl_hbr(i+1,:,:))), [-1 1])
    subplot(258)
    imagesc(corr(squeeze(med_ml_hbr(i,:,:)),squeeze(med_ml_hbr(i+1,:,:))), [-1 1])
    subplot(259)
    imagesc(corr(squeeze(med_cc_hbr(i,:,:)),squeeze(med_cc_hbr(i+1,:,:))), [-1 1])
    subplot(2,5,10)
    imagesc(corr(squeeze(med_as_hbr(i,:,:)),squeeze(med_as_hbr(i+1,:,:))), [-1 1])
    pause
    close all
    clc
    i
end


% 2 - % Differences with respect to ground truth values (by time)
% for a given density (center of the threshold range ~ 13)
diff_ge_hbo = zeros(12,14);
diff_ge_hbr = zeros(12,14);
diff_cl_hbo = zeros(12,14);
diff_cl_hbr = zeros(12,14);
diff_ml_hbo = zeros(12,14);
diff_ml_hbr = zeros(12,14);
diff_cc_hbo = zeros(12,14);
diff_cc_hbr = zeros(12,14);
diff_as_hbo = zeros(12,14);
diff_as_hbr = zeros(12,14);

dens_val = 26;

for i = 1:14
    
    temp_ge_hbo = med_ge_hbo(dens_val,:,i);
    temp_ge_hbr = med_ge_hbr(dens_val,:,i);
    temp_cl_hbo = med_cl_hbo(dens_val,:,i);
    temp_cl_hbr = med_cl_hbr(dens_val,:,i);        
    temp_ml_hbo = med_ml_hbo(dens_val,:,i);
    temp_ml_hbr = med_ml_hbr(dens_val,:,i);
    temp_cc_hbo = med_cc_hbo(dens_val,:,i);
    temp_cc_hbr = med_cc_hbr(dens_val,:,i);        
    temp_as_hbo = med_as_hbo(dens_val,:,i);
    temp_as_hbr = med_as_hbr(dens_val,:,i);
        
    for j = 1:12
        diff_ge_hbo(j,i) = (temp_ge_hbo(j)-ge_hbo_all(dens_val,i))/ge_hbo_all(dens_val,i)*100;
        diff_ge_hbr(j,i) = (temp_ge_hbr(j)-ge_hbr_all(dens_val,i))/ge_hbr_all(dens_val,i)*100;
        diff_cl_hbo(j,i) = (temp_cl_hbo(j)-cl_hbo_all(dens_val,i))/cl_hbo_all(dens_val,i)*100;
        diff_cl_hbr(j,i) = (temp_cl_hbr(j)-cl_hbr_all(dens_val,i))/cl_hbr_all(dens_val,i)*100;
        diff_ml_hbo(j,i) = (temp_ml_hbo(j)-ml_hbo_all(dens_val,i))/ml_hbo_all(dens_val,i)*100;
        diff_ml_hbr(j,i) = (temp_ml_hbr(j)-ml_hbr_all(dens_val,i))/ml_hbr_all(dens_val,i)*100;
        diff_cc_hbo(j,i) = (temp_cc_hbo(j)-cc_hbo_all(dens_val,i))/cc_hbo_all(dens_val,i)*100;
        diff_cc_hbr(j,i) = (temp_cc_hbr(j)-cc_hbr_all(dens_val,i))/cc_hbr_all(dens_val,i)*100;
        diff_as_hbo(j,i) = (temp_as_hbo(j)-as_hbo_all(dens_val,i))/as_hbo_all(dens_val,i)*100;
        diff_as_hbr(j,i) = (temp_as_hbr(j)-as_hbr_all(dens_val,i))/as_hbr_all(dens_val,i)*100;       
    end
    
end


fig1 = figure; set(fig1, 'units', 'normalized', 'outerposition', [0 0 1 1], 'Color', [1 1 1]);
subplot(231)
plot(abs(diff_ge_hbo),':','color', [0.5 0.5 0.5], 'LineWidth', 1); hold on
plot(median(abs(diff_ge_hbo),2),'color', [0.87, 0, 0.16], 'Linewidth',1.5)
ylim([0 80]);xlim([0.5 12.5]); title('Global efficiency', 'Fontsize', 20)
xlabel('Recording duration (minutes)')
ylabel('% Difference (abs)', 'Fontsize', 20)
set(gca,'Fontsize', 18); box off

%subplot(232)
%plot(abs(diff_cl_hbo),':','color', [0.5 0.5 0.5], 'LineWidth', 1); hold on
%plot(median(abs(diff_cl_hbo),2),'r', 'Linewidth',1.5)
%ylim([0 100]);xlim([0.5 12.5]); axis square; title('Modularity1')
%xlabel('Recording duration (minutes)')

subplot(232)
plot(abs(diff_ml_hbo),':','color', [0.5 0.5 0.5], 'LineWidth', 1); hold on
plot(median(abs(diff_ml_hbo),2),'color', [0.87, 0, 0.16], 'Linewidth',1.5)
ylim([0 80]);xlim([0.5 12.5]);  title('Modularity', 'Fontsize', 20)
xlabel('Recording duration (minutes)')
set(gca,'Fontsize', 18); box off

subplot(233)
plot(abs(diff_cc_hbo),':','color', [0.5 0.5 0.5], 'LineWidth', 1); hold on
plot(median(abs(diff_cc_hbo),2),'color',[0.87, 0, 0.16], 'Linewidth',1.5)
ylim([0 80]);xlim([0.5 12.5]); title('Clustering coefficient', 'Fontsize', 20)
xlabel('Recording duration (minutes)')
set(gca,'Fontsize', 18); box off
%subplot(255)
%plot(abs(diff_as_hbo),':','color', [0.5 0.5 0.5], 'LineWidth', 1); hold on
%plot(median(abs(diff_as_hbo),2),'r', 'Linewidth',1.5)
%ylim([0 100]);xlim([0.5 12.5]); axis square; title('Assortativity')
%xlabel('Recording duration (minutes)')

subplot(234)
plot(abs(diff_ge_hbr),':','color', [0.5 0.5 0.5], 'LineWidth', 1); hold on
plot(median(abs(diff_ge_hbr),2),'color', [0, 0, 0.8], 'Linewidth',1.5)
ylim([0 80]);xlim([0.5 12.5]);
xlabel('Recording duration (minutes)')
ylabel('% Difference (abs)', 'Fontsize', 20)
set(gca,'Fontsize', 18); box off

%subplot(235)
%plot(abs(diff_cl_hbr),':','color', [0.5 0.5 0.5], 'LineWidth', 1); hold on
%plot(median(abs(diff_cl_hbr),2),'b', 'Linewidth',1.5)
%ylim([0 100]);xlim([0.5 12.5]); axis square
%xlabel('Recording duration (minutes)')

subplot(235)
plot(abs(diff_ml_hbr),':','color', [0.5 0.5 0.5], 'LineWidth', 1); hold on
plot(median(abs(diff_ml_hbr),2),'color',[0, 0, 0.8], 'Linewidth',1.5)
ylim([0 80]);xlim([0.5 12.5]); 
set(gca,'Fontsize', 18); box off
xlabel('Recording duration (minutes)', 'Fontsize', 18)

subplot(236)
plot(abs(diff_cc_hbr),':','color', [0.5 0.5 0.5], 'LineWidth', 1); hold on
plot(median(abs(diff_cc_hbr),2),'color',[0, 0, 0.8], 'Linewidth',1.5)
ylim([0 80]);xlim([0.5 12.5]); 
set(gca,'Fontsize', 18); box off
xlabel('Recording duration (minutes)', 'Fontsize', 18)

%subplot(2,5,10)
%plot(abs(diff_as_hbr),':','color', [0.5 0.5 0.5], 'LineWidth', 1); hold on
%plot(median(abs(diff_as_hbr),2),'b', 'Linewidth',1.5)
%ylim([0 100]);xlim([0.5 12.5]); axis square
%xlabel('Recording duration (minutes)')











