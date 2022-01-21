function [GMSensitivity, GMmask, J_GM_wav1_norm, J_GM_wav2_norm] = DOTHUB_MakeGMSensitivityMap(J_GM_wav1,J_GM_wav2,J_vol_wav1, J_vol_wav2,SD,Thresh)

%This is Lisa's function to create a normalized GM sensitivity map, and a
%thresholded mask.

% #########################################################################
% INPUTS ##################################################################

%  J_gm =           The GM jacobian values taken from DOTHUB_Make_ToastJacobian, size
%                   nchan x nnodes GM

% J_vol =           The volume jacobian values taken from DOTHUB_Make_ToastJacobian, size
%                   nchan x nnodes volume

% SD =              The SD file for the array used in the study.  Note that
%                   the srcpos and detpos form the SD file are ignored in
%                   favour of the above defined inputs

% Thresh =          The threshold at which the sensitivity should be cut to
%                   create the binary mask


% OUTPUTS #################################################################

% GMSensitivity =   The normalized sum of the jacobian values of good channels on the GM

% GMMask =          The binary mask created by thresholding GMSensitivity
%                   at the normalized value Thresh.

% #########################################################################
% #########################################################################
% Version 0.1
% EMF, University College London, December 2019 (ho ho ho)
% #########################################################################
% #########################################################################

%First, find the maximum value of the sum of J_vol for all good channels
%(maxJvol). (BE CAREFUL AS THESE ARE ALL NEGATIVE!)

% WHAT DOES IT MEAN "GOOD CHANNELS", IS THIS AFTER PREPROCESSING??
% Keep good channels only, does considered during preprocessing
J_vol_wav1_cropped = J_vol_wav1(SD.MeasListAct(1:end/2)==1,:);
J_vol_wav2_cropped = J_vol_wav2(SD.MeasListAct(end/2+1:end)==1,:);

% Find maximum value 
MaxJvol_wav1_cropped = max(abs(J_vol_wav1_cropped),[],2);
MaxJvol_wav2_cropped = max(abs(J_vol_wav2_cropped),[],2);

% Divide J_GM by maxJvol on a channel wise basis. This should
% produce GMSensitivity
J_GM_wav1_cropped = abs(J_GM_wav1(SD.MeasListAct(1:end/2)==1,:));
J_GM_wav2_cropped = abs(J_GM_wav2(SD.MeasListAct(end/2+1:end)==1,:));
nnodes = size(J_GM_wav1_cropped,2);
J_GM_wav1_norm = J_GM_wav1_cropped./repmat(MaxJvol_wav1_cropped,1,nnodes);
J_GM_wav2_norm = J_GM_wav2_cropped./repmat(MaxJvol_wav2_cropped,1,nnodes);

% As we are now normalizing within channel, the GMSensitivity is now of a
% confusing scale.
GMSensitivity_wav1 = sum(J_GM_wav1_norm);
GMSensitivity_wav2 = sum(J_GM_wav2_norm);
GMSensitivity=(GMSensitivity_wav1+GMSensitivity_wav2)./2;

GMmask1 = any(J_GM_wav1_norm>Thresh,1);%If any channel is sensitive at both wavs, we are sensitive
GMmask2 = any(J_GM_wav2_norm>Thresh,1);
GMmask = GMmask1 & GMmask2;

