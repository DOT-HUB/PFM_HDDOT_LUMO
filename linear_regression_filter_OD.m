function filt_all = linear_regression_filter_OD (data_dod, ch_idx, lpf, sf, name, path_figures)
% The following variables should be defined in advance
% data = data (HbO or HbR) in the format time x channels
% sf = sampling frequency of the fNIRS system
% globalsignal_HbO = mean signal calculated on filtered data (dimensions = time points x 1)
% Note that global signal will be different on each participant and for HbO and HbR
% lpf = Low-pass filter

% If this doesn't work because computation takes too long try:
% 1) hpf with this method, 2) lpf (homer2), 3) GSR this method (3 step process)

% Please check data before and after filtering (e.g., time series, FC matrices, hPod)
% Please check data before and after GSR (e.g., time series, FC matrices, hPod)
% Try computing GSR in 2 ways (1 - Mean of all channels, 2 - Mean of short distance channels)
% Compute PSD of the global signal (HbO and HbR) - plot time series and PSD
% Do this for all participants

% first exclude rejected channels
data_dod = data_dod(:, ch_idx);
% 1 ? Create Legendre Polynomials (high-pass filter)
% Compute length of the dataset (samples)

n_data = size(data_dod,1);
% Compute length of the dataset (seconds)
s_data = n_data/sf;
% Calculate order of Legendre polynomials
k  = 1 + floor(s_data/150);
% Create a basis set of Legendre polynomials (L)
n = linspace(-1,1,n_data)';
L = zeros(n_data,k+1);
for i = 1:k+1
    tmp = legendre(i-1,n);
    tmp = tmp(1,:);
    L(:,i) = tmp/max(abs(tmp));
end

% 2 ? Create matrix of sines and cosines for all frequencies in the sf range
dft_matrix = dftmtx(n_data);
% Find index of the low-pass filter = frequency of interest
idx = floor((lpf/sf)*n_data);          
% Select regressors of interest
dft_matrix_lpf = dft_matrix(idx:n_data-idx+1,:);
% Select sines and cosines
sin_lpf_mtx = imag(dft_matrix_lpf);
cos_lpf_mtx = real(dft_matrix_lpf);
% Matrix (time x frequency) of sines and cosines of frequencies above lpf
lpf_mtx = [cos_lpf_mtx' sin_lpf_mtx'];

% 3 - Compute nuissance regression (1) filter and (2) filter+global signal
% (1) filter
reg_mat = [L lpf_mtx];
beta_OD = pinv(reg_mat)*data_dod;
dod_filt = data_dod - reg_mat*beta_OD;

% Reshape matrix
filt_all =  NaN(size(data_dod));
filt_all(:, ch_idx) = dod_filt;

% % (2) filter + gsr
% mHbO = mean(HbO,2);
% mHbR = mean(HbR,2);
% mHbT = mean(HbT,2);
% 
% reg_mat_HbO = [L lpf_mtx mHbO];
% reg_mat_HbR = [L lpf_mtx mHbR];
% reg_mat_HbT = [L lpf_mtx mHbT];
% 
% beta_OD = pinv(reg_mat_HbO)*HbO;
% beta_HbR = pinv(reg_mat_HbR)*HbR;
% beta_HbT = pinv(reg_mat_HbT)*HbT;
% 
% HbO_gsr = HbO - reg_mat_HbO*beta_OD;
% HbR_gsr = HbR - reg_mat_HbR*beta_HbR;
% HbT_gsr = HbT - reg_mat_HbT*beta_HbT;
% 
% % Reshape to original size and store
% conc_gsr (:,1,:) = HbO_gsr;
% conc_gsr (:,2,:) = HbR_gsr;
% conc_gsr (:,3,:) = HbT_gsr;
% 
% gsr_all =  NaN(size(HbO_orig,1),size(conc_filt,2), size(HbO_orig,2));
% for i = 1:size(conc_filt,2)
%     gsr_all(:,i, ch_idx(1:length(ch_idx)/2)) = conc_gsr(:,i,:);
% end
% 
% gsr_filt (:,1,:) = gsr_all(:,1,:);
% gsr_filt (:,2,:) = gsr_all(:,2,:);
% gsr_filt (:,3,:) = gsr_all(:,3,:);
% 
% % figure: time series + FC matrices
% fig1 = figure; set(fig1, 'units', 'normalized', 'outerposition', [0 0 1 1], 'Color', [1 1 1]);
% subplot (3,6,1:2)
% plot(HbO_filt, 'color', 'r'); hold on; plot(mHbO, 'black', 'linewidth', 1)
% box off; axis tight; ylabel('HbO', 'Fontsize', 20)
% title('Filter', 'Fontsize', 24)
% subplot (3,6,3)
% imagesc(corr(HbO_filt), [-1 1]); axis square; colorbar
% xlabel('channels')
% 
% subplot (3,6,7:8)
% plot(HbR_filt, 'color', 'b'); hold on; plot(mHbR, 'black', 'linewidth', 1)
% box off; axis tight; ylabel('HbR', 'Fontsize', 20)
% subplot (3,6,9)
% imagesc(corr(HbR_filt), [-1 1]); axis square; colorbar
% xlabel('channels')
% 
% subplot (3,6,13:14)
% plot(HbT_filt, 'color', 'g'); hold on; plot(mHbT, 'black', 'linewidth', 1)
% box off; axis tight; ylabel('HbT', 'Fontsize', 20)
% xlabel('Time (samples)', 'Fontsize', 20)
% subplot (3,6,15)
% imagesc(corr(HbT_filt), [-1 1]); axis square; colorbar
% xlabel('channels')
% 
% subplot (3,6, 4:5)
% plot(HbO_gsr, 'color', 'r'); 
% box off; axis tight
% title('Filter + GSR', 'Fontsize', 24)
% subplot (3,6,6)
% imagesc(corr(HbO_gsr), [-1 1]); axis square; colorbar; 
% xlabel('channels')
% 
% subplot (3,6, 10:11)
% plot(HbR_gsr, 'color', 'b'); 
% box off; axis tight
% subplot (3,6,12)
% imagesc(corr(HbR_gsr), [-1 1]); axis square; colorbar
% xlabel('channels')
% 
% subplot (3,6,16:17)
% plot(HbT_gsr, 'color', 'g'); 
% box off; axis tight
% xlabel('Time (samples)', 'Fontsize', 20)
% subplot (3,6,18)
% imagesc(corr(HbT_gsr), [-1 1]); axis square; colorbar
% xlabel('channels')
% colormap jet
% 
% % Save Figure
% cd(path_figures)
% saveas(fig1, [name(1:end-5) '_2filter.tiff'], 'tiffn')
% close all

end 

    
    