% ============================================================
% ======== GROUP ICA + DUAL REGRESSION (Calhoun, 2011) ========
% =============================================================
%
% In the first part of the script GROUP ICA maps/components for the selected
% PCA threshold (number of ICs) are computed
%
% Second part of the script --> FIGURES 
%
% Third part of the script --> DUAL REGRESSION
%
% Fourth part of the script --> STATS
%
% Borja Blanco 2019
% 
%
%% ------------------ First part of the script ----------------
close all; clear all; clc

% Add paths (scripts and toolboxes)
path_scripts = '/bcbl/home/home_a-f/bblanco/PROJECT_RS_SL_2017/My_scripts_RS';
addpath (genpath(path_scripts))
addpath('/bcbl/home/home_a-f/bblanco/MyMatlab/toolboxes/icasso122')
addpath('/bcbl/home/home_a-f/bblanco/MyMatlab/toolboxes/FastICA_25')

% Define experimental groups (BIL, SP, BQ)
lang_back = {'BIL'; 'SP'; 'BQ'};
ch = 46;
% Initialize group variables
oxy_std = [];
deoxy_std = []; 
sub_std3 = [];

% Concatenate data in time for temporal ICA
clc
for i = 1: length(lang_back)
    
    % Set the path to preprocessed data
    path_preprocessing =(['/bcbl/home/home_a-f/bblanco/PROJECT_RS_SL_2017/RS_SL_DATA/RS_ALL/Filter_Normal/' lang_back{i}]);
    cd(path_preprocessing)
    
    % List subjects
    sub = dir('*_preprocessed.mat');
   
    for nsub = 1:length(sub)
        
        display ([lang_back(i) sub(nsub).name(1:11)])
        
        % Load data
        data = load (sub(nsub).name);
        
        % Extract the structure --> Global Signal Regression
        %data_oxy = data.rsData.GSR_oxy;
        %data_deoxy = data.rsData.GSR_deoxy;
        
        % Extract the structure --> No global signal
        tmp_data = data.rsData.filt;
        tmp_data(:, data.rsData.CH_remove) = [];
        data_oxy = tmp_data(:, 1:ch); 
        data_deoxy = tmp_data(:, ch+1:2*ch); 

        % Normalization HbO2
        a = repmat(mean(data_oxy,1),[length(data_oxy) 1]);
        b = repmat(std(data_oxy), [length(data_oxy) 1]);
        z_oxy = (data_oxy-a)./b;
        
        % Normalization HbR
        c = repmat(mean(data_deoxy,1),[length(data_deoxy) 1]);
        d = repmat(std(data_deoxy), [length(data_deoxy) 1]);
        z_deoxy = (data_deoxy-c)./d;
       
        % Concatenate z standardized data for analysis
        oxy_std = [oxy_std; z_oxy];
        deoxy_std = [deoxy_std; z_deoxy];
        
        % Data for dual regression (data4_dualReg.mat)
        sub_std3 = cat(3, sub_std3, [z_oxy z_deoxy]);
                      
    end
    
end

hybrid_std = [oxy_std deoxy_std]';

% 1 - Dimensionality reduction (group temporal PCA)
[~ , s, latent, ~, explained] = pca(hybrid_std);

% Explained variance considering our criteria
keep_var = 0.74;

% Compute explained variance across PCA components
exp_var = cumsum(latent)./sum(latent);

% num of PCA comps with the selected explained variance
nComp = length(find(exp_var<keep_var(1)))+1;   
clc
disp(['Number of components retained = ' num2str(nComp)]);
disp(['% of variance retained =  ' num2str(exp_var(nComp)*100)]);

% Run group PCA
[coeffs, score, latent] = pca(hybrid_std,'NumComponents',nComp);
PCA_matrix = score * coeffs'; % PCA reconstructed demeaned data
PCA_matrix = bsxfun(@plus, PCA_matrix,mean(hybrid_std,2)); % plug the mean back
group_ICA = PCA_matrix; % data matrix for Icasso

% 2 - Run group ICA - ICASSO
% Force independence in the temporal domain! (channels x time)
M = 50; % number of randomizations
sR = icassoEst('randinit', group_ICA, M,'numOfIC', nComp, 'verbose', 'off', 'maxNumIterations', 500, ...
    'approach', 'defl');

% sR contains now estimates.
% Next: dissimilarity measure between them is formed
% - estimates are clustered
% - a projection for the visualization is computed

% Default similarity measure, clustering and projection:
sR = icassoExp(sR);

% Visualization & returning results
icassoShow(sR,'L',nComp,'estimate','demixing');
[iq, icaSM, w, icaTC] = icassoShow(sR,'L',nComp,'colorlimit',[.6 .7 .8 .9]);

% Order components by iq
[iq_order, idx_order] = sort(iq, 'descend');
icaSM = icaSM(:, idx_order);
icaTC = icaTC(idx_order, :);

% 3 - Store info and save data
groupICA_results.icaSM = icaSM;
groupICA_results.icaTC = icaTC;
groupICA_results.iq_order = iq_order;

% Plot correlation matrix of icaTC. 
% If they are independent, the matrix should be close to diagonal
% no correlation between components
figure
imagesc(corr(icaTC'),[-1 1])

% EXPLAINED VARIANCE BY EACH OF THE COMPONENTS (ASHBY)
% Reconstruct the signal removing one component each time
% Calculate the SSD after removing each component
% Larger values = more relative contribution of the component
% Calculate with respect to the original data
% and with respect to PCA reduced data.
ssd = zeros (2, nComp);

for i = 1: nComp
    
    % Remove one component at a time and reconstruct
    ica_temp = icaTC;
    ica_temp(i,:) = 0;
    data_rec = icaSM*ica_temp;
    
    % Compute difference
    diff_pca = group_ICA - data_rec;
    diff_orig = hybrid_std - data_rec;
    
    % Calculate sum of squares differences
    ssd(1,i) = sum(diff_pca(:).^2);
    ssd(2,i) = sum(diff_orig(:).^2);
   
end

ssd_pca = ssd(1,:)./sum(ssd(1,:))*100;
ssd_orig = ssd(2,:)./sum(ssd(2,:))*100*keep_var;

% Store info
groupICA_results.ssd_pca = ssd_pca;
groupICA_results.ssd_orig = ssd_orig;

% Correlation between oxy and deoxy spatial maps
% more negative = more consistent across oxy and deoxy
figure
oxy_Smap = icaSM(1:ch,:);
deoxy_Smap = icaSM(ch+1:2*ch,:);
corr_oxy_deoxy = corr(oxy_Smap, deoxy_Smap);
imagesc(corr_oxy_deoxy, [-1 1])
groupICA_results.corr_oxy_deoxy = diag(corr_oxy_deoxy)';

% Check similarity between spatial map extracted with ICA,
% and based on LS (regression) on the icaTC --> THEY ARE THE SAME
figure
icaSM_LS = pinv(icaTC')*group_ICA';
imagesc(icaSM_LS'-icaSM);

%% ------------------ Second part of the script -----------------
% ------------- LOAD DATA HERE IF ALREADY COMPUTED --------------
clear all; close all; clc
path_analysis ='/bcbl/home/home_a-f/bblanco/PROJECT_RS_SL_2017/RS_SL_DATA/RS_ALL/Filter_Normal/Results2019/Results/groupICA';
cd(path_analysis)
load('data4_dualReg.mat')
load('group_ICA_results.mat')

% Extract data
icaSM = groupICA_results.icaSM;
icaTC = groupICA_results.icaTC;
ssd_orig = groupICA_results.ssd_orig;
iq_order = groupICA_results.iq_order;

% Convert weights of icaSM to z-scores for plots
icaSM_std = zscore(icaSM);

% Load locations for plots
ica_locations = load('/bcbl/home/home_a-f/bblanco/PROJECT_RS_SL_2017/My_scripts_RS/ica_locations.mat');
ica_locations = ica_locations.ica_locations;

% Initialize variables for plots
ch = 46;
sf = 8.93;
N = size(icaTC,2);
freq = linspace(0, sf/2, N/2+1);
keep_var = 0.75;

% Plot components in order of Iq
for i = 1:size(icaSM,2)
    
    fig1 = figure(1); set(fig1,'units','normalized','position',[0 0 1 1], 'Color', [1 1 1])  

    % Add 10 for visualization purposes
    ica_plot = zeros (11, 8)+10;

    % Plot 1 - ICA component
    subplot(2,2,2)
    plot (icaTC(i, :)')
    title ('Time course IC  ' , 'interpreter', 'none', 'FontSize', 16)
    ylim([-10 10])
    set(gca, 'fontsize', 16)
    xlabel ('Time (samples)')
  
    % Plot 2 - single-sided amplitude spectrum.
    subplot(2,2,4)
    data  = icaTC(i, :)';
    FFT_data = abs(fft(data)); % remove the mean to improve the visualization
    bar(freq, FFT_data(1:N/2+1,:), 'blue','linewidth', 1); xlim([0, 0.1]); ylim([0 30000])
    set(gca, 'fontsize', 16)
    ylabel('Power spectrum', 'fontsize', 16)
    xlabel ('Frequency (Hz)')
    title (' FFT', 'interpreter', 'none', 'FontSize', 16)   
    
    % Component weight on each channel (oxy channels)
    temp = icaSM_std(1:ch,i);
    for loc = 1:ch
        ica_plot(ica_locations(loc,2),ica_locations(loc,3)) = temp(ica_locations(loc, 1));
    end
    
    % Plot 3 - Spatial map of the component
    subplot(2,2,1)
    plot_lim = max(abs(icaSM_std(1:ch,i)));
    %imagesc ((abs(ica_plot)>1.6).*ica_plot, [-5 5]); axis square
    imagesc (ica_plot, [-plot_lim-0.2 plot_lim+0.2]); axis square; box on; colorbar
    set(gca, 'xtick', [], 'ytick', [], 'fontsize', 16) 
    
    % Plot number of each channel
    textStrings = num2str(ica_locations(:,1));
    textStrings = strtrim(cellstr(textStrings));  % Remove any space padding
    hStrings = text(ica_locations(:,3),ica_locations(:,2),textStrings(:),...      % Plot the strings
        'HorizontalAlignment','center', 'Color', 'black', 'Fontsize', 11);
    title ('HbO_2', 'Fontsize', 16)
    
    % Component weight on each channel (deoxy channels)
    temp = icaSM_std(ch+1:2*ch,i);
    for loc = 1:ch
        ica_plot(ica_locations(loc,2),ica_locations(loc,3)) = temp(ica_locations(loc, 1));
    end
    
    % Plot 3 - spatial map of the component
    subplot(2,2,3)
    plot_lim = max(abs(icaSM_std(ch+1:2*ch,i)));
    %imagesc ((abs(ica_plot)>1.6).*ica_plot, [-5 5]); axis square
    imagesc (ica_plot, [-plot_lim-0.2 plot_lim+0.2]); axis square; box on; colorbar
    colormap hot
    set(gca, 'xtick', [], 'ytick', [], 'fontsize', 16) 
    
    % Plot number of each channel
    textStrings = num2str(ica_locations(:,1));
    textStrings = strtrim(cellstr(textStrings));  % Remove any space padding
    hStrings = text(ica_locations(:,3),ica_locations(:,2),textStrings(:),...      % Plot the strings
        'HorizontalAlignment','center', 'Color', 'black', 'Fontsize', 11);
    title ('HbR', 'interpreter', 'none', 'Fontsize', 16)
     
    % Add title
    h = suptitle(['groupICA ' num2str(i) ' - Explained variance: ' num2str(ssd_orig(i)) ' %' ]);
    set(h, 'fontsize', 16, 'interpreter', 'none')
    
    % Save figure
    screen_size = get(0, 'ScreenSize');
    origSize = get(fig1, 'Position'); % grab original on screen size
    set(fig1, 'Position', [0 0 screen_size(3) screen_size(4) ] ); % set to screen size
    set(fig1,'PaperPositionMode','auto') % set paper pos for printing
    saveas(fig1, ['Group_ICA_noGSR' num2str(i) '_hybrid_PCA' num2str(keep_var*100) '.tiff'], 'tiff')
    set(fig1,'Position', origSize) % set back to original dimensionsse(Fig1)
    close all
    
end

% Plot SPATIAL MAPS of the components only
for i = 1:size(icaSM,2)
    
    fig1 = figure(1); set(fig1,'units','normalized','position',[0 0 1 1], 'Color', [1 1 1])
    
    % add 10 for visualization purposes
    ica_plot = zeros (11, 8)+10;

    % Component weight on each channel (oxy channels)
    temp = icaSM_std(1:ch,i);
    for loc = 1:ch
        ica_plot(ica_locations(loc,2),ica_locations(loc,3)) = temp(ica_locations(loc, 1));
    end
    
    % Plot 3 - spatial map of the component
    subplot(1,2,1)
    plot_lim = max(abs(icaSM_std(1:ch,i)));
    %imagesc ((abs(ica_plot)>0.8).*ica_plot, [-plot_lim plot_lim]); axis square
    imagesc (ica_plot, [-plot_lim-0.2 plot_lim+0.2]); axis square; box on; colorbar
    colormap hot
    set(gca, 'xtick', [], 'ytick', [], 'fontsize', 16) 
    
    % Plot number of each channel
    textStrings = num2str(ica_locations(:,1));
    textStrings = strtrim(cellstr(textStrings));  % Remove any space padding
    hStrings = text(ica_locations(:,3),ica_locations(:,2),textStrings(:),...      % Plot the strings
        'HorizontalAlignment','center', 'Color', 'black', 'Fontsize', 12);
    title ('HbO_2', 'Fontsize', 16)
    
    % Component weight on each channel (deoxy channels)
    temp = icaSM_std(ch+1:2*ch,i);
    for loc = 1:ch
        ica_plot(ica_locations(loc,2),ica_locations(loc,3)) = temp(ica_locations(loc, 1));
    end
    
    % Plot 3 - spatial map of the component
    subplot(1,2,2)
    plot_lim = max(abs(icaSM_std(ch+1:2*ch,i)));
    %imagesc ((abs(ica_plot)>0.8).*ica_plot, [-plot_lim plot_lim]); axis square
    imagesc (ica_plot, [-plot_lim-0.2 plot_lim+0.2]); axis square; box on; colorbar    
    set(gca, 'xtick', [], 'ytick', [], 'fontsize', 16) 
    
    % Plot number of each channel
    textStrings = num2str(ica_locations(:,1));
    textStrings = strtrim(cellstr(textStrings));  % Remove any space padding
    hStrings = text(ica_locations(:,3),ica_locations(:,2),textStrings(:),...      % Plot the strings
        'HorizontalAlignment','center', 'Color', 'black', 'Fontsize', 12);
    title ('HbR', 'interpreter', 'none', 'Fontsize', 16)
     
    % Add title
    h = suptitle({['groupICA ' num2str(i) ' - Explained variance: ' num2str(ssd_orig(i)) ' %' ]; ['Iq = ' num2str(groupICA_results.iq_order(i))];...
        ['correlation HbO_2 - HbR = ' num2str(groupICA_results.corr_oxy_deoxy(i))]});
    set(h, 'fontsize', 18)
    
    % Save figure
    screen_size = get(0, 'ScreenSize');
    origSize = get(fig1, 'Position'); % grab original on screen size
    set(fig1, 'Position', [0 0 screen_size(3) screen_size(4) ] ); % set to screen size
    set(fig1,'PaperPositionMode','auto') % set paper pos for printing
    saveas(fig1, ['Group_ICA_components_NoGSR' num2str(i) '_hybrid_PCA' num2str(keep_var*100) '.tiff'], 'tiff')
    set(fig1,'Position', origSize) % set back to original dimensionsse(Fig1)
    
    close all
    
end

%% ------------- Third part of the script ------------------
% ----------------- DUAL REGRESSION ------------------------
clear all; close all; clc
path_analysis ='/bcbl/home/home_a-f/bblanco/PROJECT_RS_SL_2017/RS_SL_DATA/RS_ALL/Filter_Normal/Results2019/Results/groupICA';
cd(path_analysis)
load('data4_dualReg.mat')
load('group_ICA_results.mat')

% Extract data
icaSM = groupICA_results.icaSM;
icaTC = groupICA_results.icaTC;
ssd_orig = groupICA_results.ssd_orig;
iq_order = groupICA_results.iq_order;

% Convert weights of icaSM to z-scores for plots
icaSM_std = zscore(icaSM);

% Initialize variables
ch = 46;
nComp = size(icaSM_std, 2);
nsub = 99;
tp = 5000; % time points of each subject

group_sm = zeros(nComp, 2*ch, nsub);

% Remove components with low Iq or correlation
%icaSM_std(:,2) = []; 

for i = 1:nsub
    
    % The model is what determines the orientation of the matrices
    % because it has only one direction in which it is correct.
    % METHOD 1 - Two step regression
    % DATA input = preprocessed data from one subject (ch x time points)
    % MODEL input = ICA Spatial Maps from group ICA (ch x ICs)
    % BETAS output = IC associated time courses for each subject (ICs x time points)   
    Y = sub_std3(:,:,i)';
    X = icaSM_std;
    
    sub_tc = inv(X'*X)*X'*Y;
    
    % DATA input = preprocessed data from one subject (time points x ch)
    % MODEL input = IC associated time courses for each subject (time points x ICs)
    % BETAS output = IC associated spatial maps for each subject (ICs x ch)
    Y = Y';
    X = sub_tc';
    sub_sm = inv(X'*X)*X'*Y;

    % METHOD 2 - The same but based on (Erhardt, 2011) paper
    % sub_tc =  sub_std3(:,:,i) * icaSM_std' * inv((icaSM_std*icaSM_std'));
    % sub_sm = inv(sub_tc'*sub_tc)*sub_tc' * sub_std3(:,:,i);
    
    % METHOD 3 - Same outcome
    % sub_sm = linsolve(sub_tc, sub_std3(:,:,i)); 
    
    group_sm(:,:,i) = sub_sm;
    
end

% Figures spatial maps
figure
for i = 1:size(group_sm,1)
    subplot(3,5,i)
    imagesc(squeeze(group_sm(i,1:ch,:))', [-4 4]); axis square;
    colormap jet
end

figure
for i = 1:size(group_sm,1)
    subplot(3,5,i)
    imagesc(squeeze(group_sm(i,ch+1:2*ch,:))', [-4 4]); axis square;
    colormap jet
end

%% --------------- Fourth part of the script ----------------
% ---------- LOAD DATA HERE IF ALREADY COMPUTED ------------
clear all; close all; clc
path_analysis ='/bcbl/home/home_a-f/bblanco/PROJECT_RS_SL_2017/RS_SL_DATA/RS_ALL/Filter_Normal/Results2019/Results/groupICA';
cd(path_analysis)
load('dualReg_indiv_SM.mat')
load('groupICA_zmaps.mat')
ch = 46;
ica_locations = load('/bcbl/home/home_a-f/bblanco/PROJECT_RS_SL_2017/My_scripts_RS/ica_locations.mat');
ica_locations = ica_locations.ica_locations;

% Stats F-test NBS
% X model
groups = 3;
n_bil = 36; n_sp = 30; n_bq = 33;
X = zeros(n_bil+ n_sp + n_bq, groups);
X(1:n_bil,1) = 1; % BIL
X(n_bil+1:n_bil+n_sp,2) = 1; % SP
X(n_bil+n_sp +1:end,3) = 1; % BQ

% ------------------------------------------
% Exactly same outcome as below
% Prueba ANOVA
%data = squeeze(group_sm(14,1:46,:));

%data_anova = zeros(46,1);
%list_anova  = ones(99,1);
%list_anova(37:66) = 2*list_anova(37:66,1);
%list_anova(67:end) = 3*list_anova(67:end,1);
%for i = 1:46
%   [p,table, stats] = anova1(data(i,:)', list_anova);
%   clc
%   close all
%   data_anova(i) = cell2mat(table(2,5)); 
%end
% --------------------------------------------

% p = 0.05 for F(2, 96) = 3.1 uncorrected

% Contrast
c = [1 1 1];

% Number of predictors (including intercept)
p = length(c);

% Number of independent GLM's to fit (voxels)
M = size(group_sm,2);

% Number of observations (subjects)
n = size(group_sm,3);

% Stats F-test NBS
test_stat = zeros(M, size(group_sm,1));

for i= 1:size(group_sm,1)
    
    y = squeeze(group_sm(i,:,:));
    y = y';
    
    %betas = X\y;
    betas = inv(X'*X)*X'*y;
    
    %Sum of squares due to error
    sse = sum((y-X*betas).^2);
    
    %Sum of squares due to regression
    ssr = sum((X*betas-repmat(mean(y),n,1)).^2);
    
    % F test
    test_stat (:,i) = (ssr/(p-1))./(sse/(n-p));  
    
end

% figures SM for each group
bil_sm_oxy = mean(group_sm(:, 1:ch, 1:36),3);
bil_sm_deoxy = mean(group_sm(:, ch+1:2*ch, 1:36),3);
sp_sm_oxy = mean(group_sm(:, 1:ch, 37:66),3);
sp_sm_deoxy = mean(group_sm(:, ch+1:2*ch, 37:66),3);
bq_sm_oxy = mean(group_sm(:, 1:ch, 67:99),3);
bq_sm_deoxy = mean(group_sm(:, ch+1:2*ch, 67:99),3);

% Figures SM for HbO2 and HbR
% SM for the whole sample, BIL, SP, BQ and F-test
for i = 1:size(group_sm,1)
    
    fig1 = figure(1); set(fig1,'units','normalized','position',[0 0 1 1], 'Color', [1 1 1])
    
    % add 10 for visualization purposes
    ica_plotG = zeros (11, 8)+10;
    ica_plotBIL = zeros (11, 8)+10;
    ica_plotSP = zeros (11, 8)+10;
    ica_plotBQ = zeros (11, 8)+10;
    ica_plotF = zeros (11, 8)+10;
  
    % Component weight on each channel (oxy channels)
    tempG = icaSM_std(1:ch,i);
    tempBIL_oxy = bil_sm_oxy(i, :)';
    tempSP_oxy = sp_sm_oxy(i, :)';
    tempBQ_oxy = bq_sm_oxy(i, :)';
    tempF_oxy = test_stat(1:ch,i);
   
    for loc = 1:ch
        ica_plotG(ica_locations(loc,2),ica_locations(loc,3)) = tempG(ica_locations(loc, 1));
        ica_plotBIL(ica_locations(loc,2),ica_locations(loc,3)) = tempBIL_oxy(ica_locations(loc, 1));
        ica_plotSP(ica_locations(loc,2),ica_locations(loc,3)) = tempSP_oxy(ica_locations(loc, 1));
        ica_plotBQ(ica_locations(loc,2),ica_locations(loc,3)) = tempBQ_oxy(ica_locations(loc, 1));
        ica_plotF(ica_locations(loc,2),ica_locations(loc,3)) = tempF_oxy(ica_locations(loc, 1));       
    end
    
    % Plot 3 - spatial map of the component
    subplot(2,5,1)
    %plot_lim = max(abs(icaSM_std(1:ch,i)));
    %imagesc ((abs(ica_plot)>0.8).*ica_plot, [-plot_lim plot_lim]); axis square
    imagesc (ica_plotG, [-2.3 2.3]); axis square; box on; colorbar
    colormap hot
    set(gca, 'xtick', [], 'ytick', [], 'fontsize', 16) 
    
    % Plot number of each channel
    textStrings = num2str(ica_locations(:,1));
    textStrings = strtrim(cellstr(textStrings));  % Remove any space padding
    hStrings = text(ica_locations(:,3),ica_locations(:,2),textStrings(:),...      % Plot the strings
        'HorizontalAlignment','center', 'Color', 'black', 'Fontsize', 12);
    title ('HbO Group', 'Fontsize', 16)
    
    subplot(2,5,2)
    %plot_lim = max(abs(icaSM_std(1:ch,i)));
    %imagesc ((abs(ica_plot)>0.8).*ica_plot, [-plot_lim plot_lim]); axis square
    imagesc (ica_plotBIL, [-2.3 2.3]); axis square; box on; colorbar
    set(gca, 'xtick', [], 'ytick', [], 'fontsize', 16) 
    
    % Plot number of each channel
    textStrings = num2str(ica_locations(:,1));
    textStrings = strtrim(cellstr(textStrings));  % Remove any space padding
    hStrings = text(ica_locations(:,3),ica_locations(:,2),textStrings(:),...      % Plot the strings
        'HorizontalAlignment','center', 'Color', 'black', 'Fontsize', 12);
    title ('HbO BIL', 'Fontsize', 16)
        
    subplot(2,5,3)
    %plot_lim = max(abs(icaSM_std(1:ch,i)));
    %imagesc ((abs(ica_plot)>0.8).*ica_plot, [-plot_lim plot_lim]); axis square
    imagesc (ica_plotSP, [-2.3 2.3]); axis square; box on; colorbar
    set(gca, 'xtick', [], 'ytick', [], 'fontsize', 16) 
    
    % Plot number of each channel
    textStrings = num2str(ica_locations(:,1));
    textStrings = strtrim(cellstr(textStrings));  % Remove any space padding
    hStrings = text(ica_locations(:,3),ica_locations(:,2),textStrings(:),...      % Plot the strings
        'HorizontalAlignment','center', 'Color', 'black', 'Fontsize', 12);
    title ('HbO SP', 'Fontsize', 16)
        
    subplot(2,5,4)
    %plot_lim = max(abs(icaSM_std(1:ch,i)));
    %imagesc ((abs(ica_plot)>0.8).*ica_plot, [-plot_lim plot_lim]); axis square
    imagesc (ica_plotBQ, [-2.3 2.3]); axis square; box on; colorbar
    set(gca, 'xtick', [], 'ytick', [], 'fontsize', 16) 
    
    % Plot number of each channel
    textStrings = num2str(ica_locations(:,1));
    textStrings = strtrim(cellstr(textStrings));  % Remove any space padding
    hStrings = text(ica_locations(:,3),ica_locations(:,2),textStrings(:),...      % Plot the strings
        'HorizontalAlignment','center', 'Color', 'black', 'Fontsize', 12);
    title ('HbO BQ', 'Fontsize', 16)
            
    subplot(2,5,5)
    %plot_lim = max(abs(icaSM_std(1:ch,i)));
    imagesc ((abs(ica_plotF)>=3).*ica_plotF, [-3 3]); axis square; colorbar
    %imagesc (ica_plotF, [-2.3 2.3]); axis square; box on; colorbar
    colormap hot
    set(gca, 'xtick', [], 'ytick', [], 'fontsize', 16) 
    
    % Plot number of each channel
    textStrings = num2str(ica_locations(:,1));
    textStrings = strtrim(cellstr(textStrings));  % Remove any space padding
    hStrings = text(ica_locations(:,3),ica_locations(:,2),textStrings(:),...      % Plot the strings
        'HorizontalAlignment','center', 'Color', 'black', 'Fontsize', 12);
    title ('HbO F-test', 'Fontsize', 16); 
      
    % --------------------- HbR plots --------------
    % add 10 for visualization purposes
    ica_plotG = zeros (11, 8)+10;
    ica_plotBIL = zeros (11, 8)+10;
    ica_plotSP = zeros (11, 8)+10;
    ica_plotBQ = zeros (11, 8)+10;
    ica_plotF = zeros (11, 8)+10;

    % Component weight on each channel (deoxy channels)
    tempG = icaSM_std(ch+1:2*ch,i);
    tempBIL_deoxy = bil_sm_deoxy(i, :)';
    tempSP_deoxy = sp_sm_deoxy(i, :)';
    tempBQ_deoxy = bq_sm_deoxy(i, :)';
    tempF_deoxy = test_stat(ch+1:2*ch,i);

    for loc = 1:ch
        
        ica_plotG(ica_locations(loc,2),ica_locations(loc,3)) = tempG(ica_locations(loc, 1));
        ica_plotBIL(ica_locations(loc,2),ica_locations(loc,3)) = tempBIL_deoxy(ica_locations(loc, 1));
        ica_plotSP(ica_locations(loc,2),ica_locations(loc,3)) = tempSP_deoxy(ica_locations(loc, 1));
        ica_plotBQ(ica_locations(loc,2),ica_locations(loc,3)) = tempBQ_deoxy(ica_locations(loc, 1));
        ica_plotF(ica_locations(loc,2),ica_locations(loc,3)) = tempF_deoxy(ica_locations(loc, 1));
        
    end
    
    % Plot 3 - Spatial maps of the components (HbR)
    subplot(2,5,6)
    %plot_lim = max(abs(icaSM_std(ch+1:2*ch,i)));
    %imagesc ((abs(ica_plot)>0.8).*ica_plot, [-plot_lim plot_lim]); axis square
    imagesc (ica_plotG, [-2.3 2.3]); axis square; box on; colorbar
    colormap hot
    set(gca, 'xtick', [], 'ytick', [], 'fontsize', 16) 
    
    % Plot number of each channel
    textStrings = num2str(ica_locations(:,1));
    textStrings = strtrim(cellstr(textStrings));  % Remove any space padding
    hStrings = text(ica_locations(:,3),ica_locations(:,2),textStrings(:),...      % Plot the strings
        'HorizontalAlignment','center', 'Color', 'black', 'Fontsize', 12);
    title ('HbR Group', 'Fontsize', 16)
    
    subplot(2,5,7)
    %plot_lim = max(abs(icaSM_std(1:ch,i)));
    %imagesc ((abs(ica_plot)>0.8).*ica_plot, [-plot_lim plot_lim]); axis square
    imagesc (ica_plotBIL, [-2.3 2.3]); axis square; box on; colorbar
    set(gca, 'xtick', [], 'ytick', [], 'fontsize', 16) 
    
    % Plot number of each channel
    textStrings = num2str(ica_locations(:,1));
    textStrings = strtrim(cellstr(textStrings));  % Remove any space padding
    hStrings = text(ica_locations(:,3),ica_locations(:,2),textStrings(:),...      % Plot the strings
        'HorizontalAlignment','center', 'Color', 'black', 'Fontsize', 12);
    title ('HbR BIL', 'Fontsize', 16)
       
    subplot(2,5,8)
    %plot_lim = max(abs(icaSM_std(1:ch,i)));
    %imagesc ((abs(ica_plot)>0.8).*ica_plot, [-plot_lim plot_lim]); axis square
    imagesc (ica_plotSP, [-2.3 2.3]); axis square; box on; colorbar
    set(gca, 'xtick', [], 'ytick', [], 'fontsize', 16) 
    
    % Plot number of each channel
    textStrings = num2str(ica_locations(:,1));
    textStrings = strtrim(cellstr(textStrings));  % Remove any space padding
    hStrings = text(ica_locations(:,3),ica_locations(:,2),textStrings(:),...      % Plot the strings
        'HorizontalAlignment','center', 'Color', 'black', 'Fontsize', 12);
    title ('HbR SP', 'Fontsize', 16)
    
    subplot(2,5,9)
    %plot_lim = max(abs(icaSM_std(1:ch,i)));
    %imagesc ((abs(ica_plot)>0.8).*ica_plot, [-plot_lim plot_lim]); axis square
    imagesc (ica_plotBQ, [-2.3 2.3]); axis square; box on; colorbar
    set(gca, 'xtick', [], 'ytick', [], 'fontsize', 16) 
    
    % Plot number of each channel
    textStrings = num2str(ica_locations(:,1));
    textStrings = strtrim(cellstr(textStrings));  % Remove any space padding
    hStrings = text(ica_locations(:,3),ica_locations(:,2),textStrings(:),...      % Plot the strings
        'HorizontalAlignment','center', 'Color', 'black', 'Fontsize', 12);
    title ('HbR BQ', 'Fontsize', 16)
     
    subplot(2,5,10)
    %plot_lim = max(abs(icaSM_std(1:ch,i)));
    imagesc ((abs(ica_plotF)>=3).*ica_plotF, [-3 3]); axis square;colorbar
    %imagesc (ica_plotF, [-2.3 2.3]); axis square; box on; colorbar
    set(gca, 'xtick', [], 'ytick', [], 'fontsize', 16) 
    
    % Plot number of each channel
    textStrings = num2str(ica_locations(:,1));
    textStrings = strtrim(cellstr(textStrings));  % Remove any space padding
    hStrings = text(ica_locations(:,3),ica_locations(:,2),textStrings(:),...      % Plot the strings
        'HorizontalAlignment','center', 'Color', 'black', 'Fontsize', 12);
    title ('HbR F-test', 'Fontsize', 16); 
     
    % Add title
    h = suptitle(['groupICA component ' num2str(i)]);
    set(h, 'fontsize', 18)
    
    % Save figure
    screen_size = get(0, 'ScreenSize');
    origSize = get(fig1, 'Position'); % grab original on screen size
    set(fig1, 'Position', [0 0 screen_size(3) screen_size(4) ] ); % set to screen size
    set(fig1,'PaperPositionMode','auto') % set paper pos for printing
    saveas(fig1, ['Stats_groupICA_IC' num2str(i) '.tiff'], 'tiff')
    set(fig1,'Position', origSize) % set back to original dimensionsse(Fig1)
    
    close all
end



% PERMUTATION TESTING
% X model
groups = 3;
n_bil = 36; n_sp = 30; n_bq = 33;
X = zeros(n_bil+ n_sp + n_bq, groups);
X(1:36,1) = 1; % BIL
X(37:66,2) = 1; % SP
X(67:99,3) = 1; % BQ

% Contrast
c = [1 1 1];

% Number of predictors (including intercept)
p = length(c);

% Number of independent GLM's to fit (voxels)
M = size(group_sm,2);

% Number of observations (subjects)
n = size(group_sm,3);
n_perm = 10000;
test_stat = zeros(M, n_perm, size(group_sm,1));
test_threshold = zeros (2, size(group_sm,1));
for i= 1:size(group_sm,1)
    
    for j  = 1:n_perm
        
        idx_rand = randperm(99,99);
        y = squeeze(group_sm(i,:,idx_rand));
        y = y';
        
        betas = inv(X'*X)*X'*y;
        
        %Sum of squares due to error
        sse = sum((y-X*betas).^2);
        
        %Sum of squares due to regression
        ssr = sum((X*betas-repmat(mean(y),n,1)).^2);
        
        test_stat (:,j,i) = (ssr/(p-1))./(sse/(n-p));
     
    end
    test_vector = squeeze(test_stat(:,:,i));
    test_vector = sort(test_vector(:), 'ascend');
    test_idx_05 = test_vector(size(test_vector,1)*0.95);
    test_idx_01 = test_vector(size(test_vector,1)*0.99);
    test_threshold(1,i) = test_idx_05;
    test_threshold(2,i) = test_idx_01;
end


% F test to p-value
p_stats = 1 - fcdf(test_stat, 2, 96);























