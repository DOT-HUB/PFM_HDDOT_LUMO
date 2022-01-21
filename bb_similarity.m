function [group_mat_hbo, group_mat_hbr] = bb_similarity(session_fc_hbo, groundT_fc_hbo, session_fc_hbr, groundT_fc_hbr, title_lab)

% Compute groundT network 
% For debugging
% groundT_fc_hbo = group_visC_hbo;
% session_fc_hbo = visC_hbo;
% -------------------------

% Create different groundT for each session
tp = size(session_fc_hbo,1);
ini_ses = 1:tp:size(groundT_fc_hbo,1); 
end_ses = tp:tp:size(groundT_fc_hbo,1);

% Create upper triangular mask
mask = triu(true(size(groundT_fc_hbo,2),size(groundT_fc_hbo,2)),1); 

groundT_ses = [];
groundT_fc = [];
groundT_hbo = [];
for i = 1:size(session_fc_hbo,3)
    
    % Concatenate all sessions except the actual session - JU
    temp = groundT_fc_hbo;
    temp(ini_ses(i):end_ses(i),:) = [];
    groundT_ses(:,:,i) = temp;
    % Compute FC matrix of ground T data for each session
    groundT_fc (:,:,i) = atanh(corr(groundT_ses(:,:,i)));
    % Create upper triangular mask
    temp_fc = groundT_fc(:,:,i);
    groundT_hbo(:,i) = temp_fc(mask);

end

% Compute each session fc
% session_fc_hbo = visC_hbo; % for debugging
session_indiv_hbo = [];
for i = 1:size(session_fc_hbo,3)   
    temp = atanh(corr(session_fc_hbo(:,:,i)));
    session_indiv_hbo(:,i) = temp(mask);
end

% Compute similarity between sessions and with groundT
group_mat_hbo = zeros(size(session_fc_hbo,3)+1, size(session_fc_hbo,3)+1);
group_mat_hbo(1:size(session_fc_hbo,3), 1:size(session_fc_hbo,3)) = corr(session_indiv_hbo);
group_mat_hbo(size(session_fc_hbo,3)+1, size(session_fc_hbo,3)+1) = 1;

for i = 1:size(session_fc_hbo,3)
    r_ses = corr([session_indiv_hbo(:,i), groundT_hbo(:,i)]);
    group_mat_hbo(i,size(session_fc_hbo,3)+1) = r_ses(1,2);
    group_mat_hbo(size(session_fc_hbo,3)+1, i) = r_ses(1,2);
end

% Same procedure for HBR
groundT_ses = [];
groundT_fc = [];
groundT_hbr = [];
for i = 1:size(session_fc_hbo,3)
    
    % Concatenate all sessions except the actual session - JU
    temp = groundT_fc_hbr;
    temp(ini_ses(i):end_ses(i),:) = [];
    groundT_ses(:,:,i) = temp;
    % Compute FC matrix of ground T data for each session
    groundT_fc (:,:,i) = atanh(corr(groundT_ses(:,:,i)));
    % Create upper triangular mask
    temp_fc = groundT_fc(:,:,i);
    groundT_hbr(:,i) = temp_fc(mask);

end

% Compute each session fc
session_indiv_hbr = [];
for i = 1:size(session_fc_hbr,3)    
    temp = atanh(corr(session_fc_hbr(:,:,i)));
    session_indiv_hbr(:,i) = temp(mask);
end

% Compute similarity between sessions and with groundT
group_mat_hbr = zeros(size(session_fc_hbr,3)+1, size(session_fc_hbr,3)+1);
group_mat_hbr(1:size(session_fc_hbr,3), 1:size(session_fc_hbr,3)) = corr(session_indiv_hbr);
group_mat_hbr(size(session_fc_hbr,3)+1, size(session_fc_hbr,3)+1) = 1;

for i = 1:size(session_fc_hbr,3)
    r_ses = corr([session_indiv_hbr(:,i), groundT_hbr(:,i)]);
    group_mat_hbr(i,size(session_fc_hbr,3)+1) = r_ses(1,2);
    group_mat_hbr(size(session_fc_hbr,3)+1, i) = r_ses(1,2);
end

% Plot outcomes
cmap = load('/Users/borjablanco/Documents/GitHub/Hyper24RSFC_scripts/scripts/BB_scripts/colormats/lajolla');
fig1 = figure; set(fig1, 'units', 'normalized', 'outerposition', [0 0 1 1], 'Color', [1 1 1]);

cmap = cmap.lajolla;

% HbO
subplot(121)
imagesc(group_mat_hbo, [0.4, 1])
axis square
ses_ticks = 1:1:size(session_fc_hbo,3)+1;
ses_lab = {'S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S7', ...
    'S8', 'S9', 'S10', 'S11', 'S12', 'S13', 'S14', 'GroundT'};
set(gca,'Fontsize', 28, 'YTick', ses_ticks, 'YTickLabel', ses_lab, 'XTick', ses_ticks,...
    'XTickLabel', ses_lab, 'XTickLabelRotation', 90)
title([title_lab], 'FontSize', 36, 'FontWeight', 'normal')

% HbR
subplot(122)
imagesc(group_mat_hbr, [0.4, 1])
colormap(cmap)
title([title_lab])
axis square

ses_ticks = 1:1:size(session_fc_hbo,3)+1;
ses_lab = {'S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S7', ...
    'S8', 'S9', 'S10', 'S11', 'S12', 'S13', 'S14', 'GroundT'};
set(gca,'Fontsize', 28, 'YTick', ses_ticks, 'YTickLabel', ses_lab, 'XTick', ses_ticks,...
    'XTickLabel', ses_lab, 'XTickLabelRotation', 90)
c=colorbar;
set(c, 'Fontsize', 40)
x=get(c,'Position'); x(3)=0.025; x(1)=0.92; set(c,'Position',x)
title([title_lab], 'FontSize', 36, 'FontWeight', 'normal')

saveas(fig1, ['similarity_lajolla' title_lab], 'tiff')













