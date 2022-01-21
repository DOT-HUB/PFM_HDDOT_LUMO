
function bb_similarity_plot (net_mat_hbo, net_mat_hbr, nparc, nsub, title_lab)

% Load colormap and create figure
load('/Users/borjablanco/Documents/GitHub/Hyper24RSFC_scripts/scripts/BB_scripts/colormats/roma')
fig1 = figure; set(fig1, 'units', 'normalized', 'outerposition', [0 0 1 1], 'Color', [1 1 1]);

% Plot Fisher's transform of FC matrix
% HbO
subplot(121)
imagesc(corr(net_mat_hbo), [-1 1]); hold on

for i = 1:nsub
    line([0 size(net_mat_hbo,2)], [nparc+0.5+nparc*(i-1), nparc+0.5+nparc*(i-1)],'color', 'k', 'LineWidth', 1);
    line([nparc+0.5+(nparc*(i-1)), nparc+0.5+(nparc*(i-1))], [0 size(net_mat_hbo,2)],'color', 'k', 'LineWidth', 1);
end
ses_ticks = nparc/2 + 0.5:nparc:size(net_mat_hbo,2);
ses_lab = {'S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S7', ...
    'S8', 'S9', 'S10', 'S11', 'S12', 'S13', 'S14'};
set(gca,'Fontsize', 24, 'YTick', ses_ticks, 'YTickLabel', ses_lab, 'XTick', ses_ticks,...
    'XTickLabel', ses_lab, 'XTickLabelRotation', 90)
axis square;
title([title_lab ' HbO'])

subplot(122)
imagesc(corr(net_mat_hbr), [-1 1]); hold on

for i = 1:nsub
    line([0 size(net_mat_hbo,2)], [nparc+0.5+nparc*(i-1), nparc+0.5+nparc*(i-1)],'color', 'k', 'LineWidth', 1);
    line([nparc+0.5+(nparc*(i-1)), nparc+0.5+(nparc*(i-1))], [0 size(net_mat_hbo,2)],'color', 'k', 'LineWidth', 1);
end
ses_ticks = nparc/2 + 0.5:nparc:size(net_mat_hbo,2);
ses_lab = {'S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S7', ...
    'S8', 'S9', 'S10', 'S11', 'S12', 'S13', 'S14'};
set(gca,'Fontsize', 24, 'YTick', ses_ticks, 'YTickLabel', ses_lab, 'XTick', ses_ticks,...
    'XTickLabel', ses_lab, 'XTickLabelRotation', 90)
axis square; 
title([title_lab ' HbR'])
colormap(flipud(roma)); c=colorbar;
x=get(c,'Position'); x(3)=0.025; x(1)=0.92; set(c,'Position',x)

saveas(fig1, ['similarity_allParcels' title_lab], 'tiff')
%close

% subplot(133)
% neg_mat = atanh(corr([net_mat_hbo net_mat_hbr]));
% neg_mat = neg_mat(size(net_mat_hbo,2)+1:end, 1:size(net_mat_hbo,2));
% 
% imagesc(neg_mat, [-0.75 0.75])
% for i = 1:nsub
%     line([0 size(net_mat_hbo,2)], [size1+0.5+size1*(i-1), size1+0.5+size1*(i-1)],'color', 'k', 'LineWidth', 1);
%     line([size1+0.5+(size1*(i-1)), size1+0.5+(size1*(i-1))], [0 size(net_mat_hbo,2)],'color', 'k', 'LineWidth', 1);
% end
% set(gca,'Fontsize', 24, 'YTick', ses_ticks, 'YTickLabel', ses_lab, 'XTick', ses_ticks,...
%     'XTickLabel', ses_lab, 'XTickLabelRotation', 90)
% axis square



