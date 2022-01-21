%HyperRSFCnirs2preprocessed

% Pipeline to preprocesses Rob's hyperscanning .nirs files for RSFC
% analyses. This script assumes that the original .nirs files have been
% cropped to include the rest period only.

%-----------------------Main Steps------------------------------------
%   1. Run enPrune on the raw intensity data;
%   2. Convert to DOD
%   3. DOD to conc (HbO, HbR, HbT)
%   4A. Filter to RSFC range (linear regression filter)
%   4B. Short seperation regression (Method 4)
%   6. Convert to concentration (HbO, HbR)
%   7. Perform global signal regression on concentration data (lpf = 0.08)
%   8. Convert the concentration data back to dod for reconstruction
%   9. Make a .prepro file (for use in DOTHUB-toolbox reconstruction)
%   10. Save all of the above data in a '_REST_preprocessed.mat' file
clear; close all; clc

path_preprocessing = '/Volumes/CAM_data/neoLAB/preprocessed';
path_nirs = '/Volumes/CAM_data/neoLAB/Rest_nirs_files';
path_sd3d = '/Volumes/CAM_data/neolab/SD3D_files';
path_figures = '/Volumes/CAM_data/neoLAB/indiv_figures';

path_scripts = '/Users/borjablanco/Documents/GitHub/Hyper24RSFC_scripts/scripts/BB_scripts';
addpath(genpath(path_scripts));

path_homer2 = '/Users/borjablanco/Documents/MATLAB/toolboxes/homer2';
addpath(genpath(path_homer2));

% path_wavelets = '/Users/borjablanco/Documents/MATLAB/toolboxes/BrainWavelet_v1.0_Matrix_OSX_Mavericks';
% addpath(genpath(path_wavelets));
% setup

path_DOTHUB = '/Users/borjablanco/Documents/GitHub/DOT-HUB_toolbox';
addpath(genpath(path_DOTHUB));

% ---- Wavelet parameters ----
%thresh = 0.02;
%nscale = 'extreme';
%chsear = 'moderate';
% ----------------------------

% Dependencies:
%   1. Homer toolbox
%   2. DOT-HUB_toolbox
%   3. julie GSR filters

% Updates:
% Dec 2020 JU added SSR step after regular filtering. Culver used this in infants
% after filtering
% Dec 19 2020, JU changed so that final step converted concentration back
% to dod after GSR filtering. Found that could not filter dod with GSR
% filter. See the figures comparing the three 'Hyper24_n_GSRatDifferentSteps

cd(path_sd3d);
sd3d_files = dir('*.SD3D');

cd(path_nirs);
files = dir('*.nirs');

tic;
for nsub = 1:length(files)
    
    % Load data HOMER2 format .nirs
    cd(path_nirs);
    data = load(files(nsub).name,'-mat');
    data.name = files(nsub).name;
    
    % Update SD3D file
    cd(path_sd3d)
    load(sd3d_files(nsub).name,'-mat');
    data.SD3D = SD3D;
    
    % Add relevant info to workspace and new data structure
    data.nchannels = size(data.d,2);
    data.ppf = [6.48 5.4]; % Adapted for an adult
    data.sf = abs(1/(data.t(2)-data.t(1))); % Sampling rate
    data.time = (0:size(data.d,1)-1)/data.sf; % time vector (seconds)
    data.total_length = (data.time(end)/60);
    data.lpf = 0.08;

    display([num2str(data.total_length) ' minutes'])
    
    % Take out any negative or 0 intensity values
    data.d(data.d<=0) = min(data.d(data.d>0));
    % -----------------------------------------------------------------------
    % STEP 1: Prune channels using homer2 enPruneChannels
    SNRthresh = 12.0;
    SDrange = [0 60];
    tInc = ones((size(data.d,1)),1);
    dRange = [5e-5 2.5];
    data.SD3D = enPruneChannels(data.d,data.SD3D,tInc,dRange,SNRthresh,SDrange,true);
    data.ch_idx = find(data.SD3D.MeasListAct==1);          
    %-------------------------------------------------------------------------
    % STEP 2: Convert intensity to OD values
    data.dod = hmrIntensity2OD(data.d);    
    %--------------------------------------------------------------------------
    % STEP 3: Convert OD data to concentration data.
    data.conc = hmrOD2Conc(data.dod, data.SD3D, data.ppf);
    %--------------------------------------------------------------------------
    % PREPROCESSING A
    % STEP 4A: Filter the OD data to RSFC range using linear regression filter
    %[data.conc_filt, data.conc_filtgsr] = linear_regression_filter(data.conc, data.ch_idx, data.lpf, data.sf, data.name, path_figures);   
    % --------------------------------------------------------------------    
    % PREPROCESSING B
    % STEP 4B1: Short separation regression on conc data (dc data)
    % Apply short separation regression regressing out short channels 
    % using the average of the short separation channels closer ( < rhoSD_ssThresh) 
    % to the location of the source and detector of a long Src-Det separation channel
    data.flagSSmethod = 2;
    data.rhoSD_ssThresh = 10;
    data.conc_ssr = DOTHUB_hmrSSRegressionByChannel(data.conc, data.SD3D, data.rhoSD_ssThresh, data.flagSSmethod);
    
    % STEP 4B2: Filter using the same filter as A
    [data.conc_filtssr, ~] = linear_regression_filter(data.conc_ssr, data.ch_idx, data.lpf, data.sf, data.name, path_figures);
    %--------------------------------------------------------------------------       
    % PREPROCESSING A & B comparison plots
    %bb_prep_plots(data, path_figures)
    % -------------------------------------------------------------------------
    % STEP 5A: Convert Conc to OD for reconstruction,, use hmrConc2OD
    %data.dod_filt = hmrConc2OD(data.conc_filt, data.SD3D, data.ppf);
    %data.dod_filtgsr = hmrConc2OD(data.conc_filtgsr, data.SD3D, data.ppf);
    data.dod_filtssr = hmrConc2OD(data.conc_filtssr, data.SD3D, data.ppf);
    %--------------------------------------------------------------------------
    % STEP 6: Downsample to 1/4th the original frequency (5.26/4 = 1.315 Hz)     
    %data.dod_filt_dwn = downsample(data.dod_filt,4);
    %data.dod_filtgsr_dwn = downsample(data.dod_filtgsr,4);
    data.dod_filtssr_dwn = downsample(data.dod_filtssr,4);
    data.time_dwn = downsample(data.time,4);
    % ------------------------------------------------------------------------
    % STEP 7: make a prepro file for use in the DOT-HUB_toolbox
    ds = datestr(now,'yyyymmDDHHMMSS');
    nirsFileName=([path_nirs, '/' files(nsub).name]); 
    [pathstr,name,~] = fileparts(nirsFileName);
    %preproFileNameA = [nirsFileName(1:end-5) '_filt.prepro']; % filter
    %preproFileNameB = [nirsFileName(1:end-5) 'gsr.prepro']; % filter + gsr
    preproFileNameC = [nirsFileName(1:end-5) '_ssr_avg.prepro']; % filter + ssr

    logData(1,:) = {'Created on: '; ds};
    logData(2,:) = {'Derived from data: ', nirsFileName};
    logData(3,:) = {'Pre-processed using:', mfilename('fullpath')};
    
    % write and save .prepro file. This is used during reconstruction in the DOTHUB
    % toolbox. will automatically save the file to where the nirs file came from.
    %tDOD = size(data.d,1); % this should be: The time vector in seconds corresponding to the first dimension of dod.
    %[preproA, preproFileNameA] = DOTHUB_writePREPRO(preproFileNameA,logData,data.dod_filt_dwn,data.time_dwn,data.SD3D);
    %[preproB, preproFileNameB] = DOTHUB_writePREPRO(preproFileNameB,logData,data.dod_filtgsr_dwn,data.time_dwn,data.SD3D);
    [preproC, preproFileNameC] = DOTHUB_writePREPRO(preproFileNameC,logData,data.dod_filtssr_dwn,data.time_dwn,data.SD3D);
    %-------------------------------------------------------------------------
    % STEP 8: Save as snirf
       
end
toc



