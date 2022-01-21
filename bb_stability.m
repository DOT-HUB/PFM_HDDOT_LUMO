function bb_stability (group_net_hbo, groundT_fc_hbo, group_net_hbr, groundT_fc_hbr, title_lab)

% For debugging ---------------------
% group_net_hbo = visC_hbo;
% group_net_hbr = visC_hbr;

% groundT_fc_hbo = group_visC_hbo;
% groundT_fc_hbr = group_visC_hbr;
% -----------------------------------

% Compute ground truth FC matrix based on concatenated data
% For parcels belonging to each specific network
% HbO
% Fisher transform for normality
% Create different groundT for each session
tp = 947;
ini_ses = 1:tp:size(groundT_fc_hbo,1); 
end_ses = tp:tp:size(groundT_fc_hbo,1);

% Create upper triangular mask
mask = triu(true(size(groundT_fc_hbo,2),size(groundT_fc_hbo,2)),1); 

groundT_ses = [];
groundT_fc = [];
groundT_precFC_hbo = [];
for i = 1:size(group_net_hbo,3)
    
    % Concatenate all sessions except the actual session - JU
    temp = groundT_fc_hbo;
    temp(ini_ses(i):end_ses(i),:) = [];
    groundT_ses(:,:,i) = temp;
    % Compute FC matrix of ground T data for each session
    groundT_fc (:,:,i) = atanh(corr(groundT_ses(:,:,i)));
    % Create upper triangular mask
    temp_fc = groundT_fc(:,:,i);
    groundT_precFC_hbo(:,i) = temp_fc(mask);

end


% Same procedure for HBR
groundT_ses = [];
groundT_fc = [];
groundT_precFC_hbr = [];
for i = 1:size(group_net_hbo,3)
    
    % Concatenate all sessions except the actual session - JU
    temp = groundT_fc_hbr;
    temp(ini_ses(i):end_ses(i),:) = [];
    groundT_ses(:,:,i) = temp;
    % Compute FC matrix of ground T data for each session
    groundT_fc (:,:,i) = atanh(corr(groundT_ses(:,:,i)));
    % Create upper triangular mask
    temp_fc = groundT_fc(:,:,i);
    groundT_precFC_hbr(:,i) = temp_fc(mask);

end

% Define parameters for sliding window
sf = 1.315; % downsampled from 5.26Hz/4 = 1.315
data_size = 936; % samples, to have integer outputs below
% define time windows (around 60 seconds time windows)
% 1,2,3,4,5,6,7,8,9,10,11,12 minutes
tw = 1:12;
twp = data_size/length(tw);
tw_vect = tw.*twp;

% Initialize variables
prec_hbo = [];
prec_hbr = [];

for nsub = 1:size(group_net_hbo,3)
    nsub
    % Reduce data to 936 samples
    z_hbo = group_net_hbo(1:twp*length(tw),:,nsub)';
    z_hbr = group_net_hbr(1:twp*length(tw),:,nsub)';

    % Analysis: 1 to 12 minutes with 30 seconds overlap
    [prec_hbo(nsub, :), prec_hbr(nsub,:), tw_order] = bb_fcTimeWindow(z_hbo, z_hbr, tw_vect, mask, groundT_precFC_hbo(:,nsub), groundT_precFC_hbr(:, nsub));

end

% Compute median similarity in each individual for each time window
med_hbo = [];
med_hbr = [];
ntw = [];
for i = 1:length(tw)
    tw_idx = find(tw_order==i);
    med_hbo(:,i) = median(prec_hbo(:, tw_idx),2);
    med_hbr(:,i) = median(prec_hbr(:, tw_idx),2);
    ntw(i) = length(tw_idx);
end

% PLOTS
fig1 = figure; set(fig1, 'units', 'normalized', 'outerposition', [0 0 1 1], 'Color', [1 1 1]);
cmap = load('/Users/borjablanco/Documents/GitHub/Hyper24RSFC_scripts/scripts/BB_scripts/colormats/lajolla');
cmap = cmap.lajolla;
subplot(2,5,1:3)
x_ticks_tw = [ntw(1)/2 cumsum(ntw)+ [ntw(2:end)/2 0]] +0.5;
x_ticks_ds = x_ticks_tw(1:end-2);
x_label_tw = {'1', '2', '3', '4','5', '6', '7', '8', '9', '10', '11'};
y_label_sub = {'1', '2', '3', '4','5', '6', '7', '8', '9', '10', '11',  '12', '13', '14'};
imagesc(prec_hbo, [0 1]);

hold on
idx_ntw = ntw(1);
for i =1:length(ntw)-1 
    line([idx_ntw+0.5 idx_ntw+0.5], [0 14.5], 'color', 'k', 'LineWidth', 0.5)
    idx_ntw = idx_ntw + ntw(i+1);
end

for i =1:14
    line([0 144.5], [i+0.5 i+0.5], 'color', 'k', 'LineWidth', 0.5)
end
%title('HbO'); 
ylabel('Session number', 'Fontsize', 18)
set(gca, 'Fontsize',18,'Xtick',x_ticks_ds, 'XTickLabel', x_label_tw, 'XTickLabelRotation',0,'Ytick', 1:1:14, 'YTickLabel', y_label_sub);

subplot(2,5,6:8)
imagesc(prec_hbr, [0 1]);
hold on
idx_ntw = ntw(1);
for i =1:length(ntw)-1   
    line([idx_ntw+0.5 idx_ntw+0.5], [0 14.5] ,'color', 'k', 'LineWidth', 0.5)
    idx_ntw = idx_ntw + ntw(i+1);
end
for i =1:14
    line([0 144.5], [i+0.5 i+0.5], 'color', 'k', 'LineWidth', 0.5)
end
%title('HbR', 'Fontsize', 18); 
ylabel('Session number', 'Fontsize', 18)
xlabel('Time windows (minutes)', 'Fontsize', 18)
set(gca, 'Fontsize',18,'Xtick',x_ticks_ds, 'XTickLabel', x_label_tw, 'XTickLabelRotation',0,'Ytick', 1:1:14, 'YTickLabel', y_label_sub);
colormap(cmap)

% Compute and plot percentage change in right y axis

temp_hbo = median(med_hbo);
temp_hbr = median(med_hbr);

pdiff_hbo = [];
pdiff_hbr = [];

for i = 1:length(temp_hbo)-1
    pdiff_hbo(i) = (temp_hbo(i+1) - temp_hbo(i))/temp_hbo(i)*100;
    pdiff_hbr(i) = (temp_hbr(i+1) - temp_hbr(i))/temp_hbr(i)*100;

end

subplot(2,5,4:5)
yyaxis right
plot(med_hbo','color',  [0.87, 0, 0.16], 'LineStyle', ':', 'LineWidth', 1, 'Marker', 'none'); hold on
plot(median(med_hbr), 'color', [0.87, 0, 0.16], 'LineWidth', 2, 'Marker', 'none');set(gca, 'Fontsize', 18)
ylabel('Reliability (R^2)', 'Fontsize', 18)
ylim([0 1])
ax = gca;
ax.YColor = [0.87, 0, 0.16];
yyaxis left
area(1.5:1:11.5,pdiff_hbo, 'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0.5 .5 .5],'LineWidth',1);
ylabel('% Difference', 'Fontsize', 18)
xlabel('Recording duration (minutes)', 'Fontsize', 18)
ylim([0 55])
ax = gca;
ax.YColor = 'black';
xlim([0.5 12.5]);
box off
subplot(2,5,9:10)
yyaxis right
plot(med_hbr','color',  [0, 0, 0.8], 'LineStyle', ':', 'LineWidth', 1, 'Marker', 'none'); hold on
plot(median(med_hbr),'color', [0, 0, 0.8], 'LineWidth', 2, 'Marker', 'none');set(gca, 'Fontsize', 18)
ylabel('Reliability (R^2)', 'Fontsize', 18)
ylim([0 1])
ax = gca;
ax.YColor = [0, 0, 0.8];
yyaxis left
area(1.5:1:11.5,pdiff_hbr, 'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0.5 .5 .5],'LineWidth',1);
ylabel('% Difference', 'Fontsize', 18)
xlabel('Recording duration (minutes)', 'Fontsize', 18)
ylim([0 55])
xlim([0.5 12.5]);
ax = gca;
ax.YColor = 'black';
box off

saveas(fig1, ['stability_' title_lab], 'tiff')

