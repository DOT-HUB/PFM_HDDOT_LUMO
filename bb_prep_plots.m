function bb_prep_plots(data, path_figures)

% Create variables for plots
ch_idx = data.ch_idx(1:size(data.ch_idx,1)/2);

HbO_filt = squeeze(data.conc_filt(:,1,ch_idx));
HbR_filt = squeeze(data.conc_filt(:,2,ch_idx));
% (2) filter + gsr
mHbO = mean(HbO_filt,2);
mHbR = mean(HbR_filt,2);

HbO_gsr = squeeze(data.conc_filtgsr(:,1,ch_idx));
HbR_gsr = squeeze(data.conc_filtgsr(:,2,ch_idx));

HbO_ssr = squeeze(data.conc_filtssr(:,1,ch_idx));
HbR_ssr = squeeze(data.conc_filtssr(:,2,ch_idx));

% figure: time series + FC matrices
fig1 = figure; set(fig1, 'units', 'normalized', 'outerposition', [0 0 1 1], 'Color', [1 1 1]);

subplot (3,6,1:2)
plot(data.time, HbO_filt, 'color', 'r'); hold on; plot(data.time, mHbO, 'black', 'linewidth', 1)
title ('HbO', 'Fontsize', 16)
box off; axis tight; ylabel({'Filter' ; 'Conc. (\muM)'}, 'Fontsize', 14)
ylim([-max(abs(HbO_filt(:))) max(abs(HbO_filt(:)))])

subplot (3,6,3)
imagesc(corr(HbO_filt), [-1 1]); axis square; colorbar
xlabel('Channels', 'Fontsize', 12)

subplot (3,6,4:5)
plot(data.time, HbR_filt, 'color', 'b'); hold on; plot(data.time, mHbR, 'black', 'linewidth', 1)
title ('HbR', 'Fontsize', 16)
box off; axis tight;
ylim([-max(abs(HbR_filt(:))) max(abs(HbR_filt(:)))])

subplot (3,6,6)
imagesc(corr(HbR_filt), [-1 1]); axis square; colorbar
xlabel('Channels', 'Fontsize', 12)

subplot (3,6,7:8)
plot(data.time, HbO_gsr, 'color', 'r'); 
box off; axis tight; ylabel({'Filter-gsr' ; 'Conc. (\muM)'}, 'Fontsize', 14)
ylim([-max(abs(HbO_gsr(:))) max(abs(HbO_gsr(:)))])

subplot (3,6,9)
imagesc(corr(HbO_gsr), [-1 1]); axis square; colorbar
xlabel('Channels', 'Fontsize', 12)

subplot (3,6, 10:11)
plot(data.time, HbR_gsr, 'color', 'b'); 
box off; axis tight
ylim([-max(abs(HbR_gsr(:))) max(abs(HbR_gsr(:)))])

subplot (3,6,12)
imagesc(corr(HbR_gsr), [-1 1]); axis square; colorbar; 
xlabel('Channels', 'Fontsize', 12)

subplot (3,6, 13:14)
plot(data.time, HbO_ssr, 'color', 'r');
box off; axis tight; ylabel({'Filter-ssr' ; 'Conc. (\muM)'}, 'Fontsize', 14)
xlabel('Time (seconds)', 'Fontsize', 12)
ylim([-max(abs(HbO_ssr(:))) max(abs(HbO_ssr(:)))])

subplot (3,6,15)
imagesc(corr(HbO_ssr), [-1 1]); axis square; colorbar
xlabel('Channels', 'Fontsize', 12)

subplot (3,6,16:17)
plot(data.time, HbR_ssr, 'color', 'b'); 
box off; axis tight
xlabel('Time (seconds)', 'Fontsize', 12)
ylim([-max(abs(HbR_ssr(:))) max(abs(HbR_ssr(:)))])

subplot (3,6,18)
imagesc(corr(HbR_ssr), [-1 1]); axis square; colorbar
xlabel('Channels', 'Fontsize', 12)
colormap jet

% Save Figure
cd(path_figures)
saveas(fig1, [data.name(1:end-5) '_filter_comp.tiff'], 'tif')
close all


