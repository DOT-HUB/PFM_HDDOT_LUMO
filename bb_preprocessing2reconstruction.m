
% Registration script
clear; close all; clc

path_toast = '/Users/borjablanco/Documents/MATLAB/toolboxes/toast';
addpath(genpath(path_toast));

path_toast2 = '/Users/borjablanco/Documents/MATLAB/toolboxes/toastpp-2.0.2';
addpath(genpath(path_toast2));

path_DOTHUB = '/Users/borjablanco/Documents/GitHub/DOT-HUB_toolbox';
addpath(genpath(path_DOTHUB));

path_homer2 = '/Users/borjablanco/Documents/MATLAB/toolboxes/homer2';
addpath(genpath(path_homer2));

path_sd3d = '/Volumes/CAM_data/neolab/SD3D_files';
path_prepro = '/Volumes/CAM_data/neolab/Rest_nirs_files';
mesh_path = '/Volumes/CAM_data/neolab/Rob_Mesh/Robmesh_bb.mshs';

%load(mesh_path, '-mat')

% fileName = '/Volumes/CAM_data/neolab/Rob_Mesh/Robmesh_bb.mshs';
% rob_mesh.fileName = fileName;
% rob_mesh.landmarks = landmarks;
% rob_mesh.gmSurfaceMesh = gmSurfaceMesh;
% rob_mesh.headVolumeMesh = headVolumeMesh;
% rob_mesh.scalpSurfaceMesh = scalpSurfaceMesh;
% rob_mesh.vol2gm = vol2gm;
% 
% save('Robmesh_bb.mshs', '-struct', 'rob_mesh')


% ============================ START HERE
% ====================== IF COMPUTING JACOBIAN
cd(path_sd3d);
files_sd3d = dir('*.SD3D');

for nsub = 1:length(files_sd3d)  
    
    SD3DFileName = files_sd3d(nsub).name;
    % Register chosen mesh to subject SD3D and create rmap
    [rmap, rmapFileName] = DOTHUB_meshRegistration(SD3DFileName, mesh_path);
    figure;
    DOTHUB_plotRMAP(rmap)
    
end

% Compute Jacobian
cd(path_sd3d);
files_rmap = dir('*.rmap');

% This has been computed on Bodhi
for nsub = 1:length(files_rmap)
    
    tic;
    % Calculate Jacobian
    basis = [30 30 30];
    rmapFileName = files_rmap(nsub).name;
    [jac, jacFileName] = DOTHUB_makeToastJacobian(rmapFileName,basis);    
    toc
    
end

% ================================ START HERE
% ====================== IF JACOBIAN ALREADY COMPUTED
path_DOTHUB = '/Users/borjablanco/Documents/GitHub/DOT-HUB_toolbox';
addpath(genpath(path_DOTHUB));
path_toast = '/Users/borjablanco/Documents/MATLAB/toolboxes/toast';
addpath(genpath(path_toast));
path_toast2 = '/Users/borjablanco/Documents/MATLAB/toolboxes/toastpp-2.0.2';
addpath(genpath(path_toast2));

path_sd3d = '/Volumes/CAM_data/neolab/SD3D_files';
path_prepro = '/Volumes/CAM_data/neolab/Rest_nirs_files';

% Invert Jacobian
cd(path_sd3d);
files_jac = dir('*.jac');
cd(path_prepro);
files_prepro = dir ('*ssr.prepro');

% In case it needs renaming
% for nsub = 1:length(files_jac)
%     jacFileName = [path_sd3d filesep files_jac(nsub).name];
%     load(jacFileName, 'fileName', '-mat');
%     fileName = jacFileName
%     save(jacFileName, 'fileName', '-append')
% end

% Invert Jacobian
for nsub = 1:length(files_jac)
    tic;
    jacFileName = [path_sd3d filesep files_jac(nsub).name];
    preproFileName = [path_prepro filesep files_prepro(nsub).name];
    %Note that you can either separately calculate the inverse, or just run
    %DOTHUB_reconstruction, which will then call the inversion itself
    [invjac, invjacFileName] = DOTHUB_invertJacobian(jacFileName,preproFileName,'saveFlag',true);
    toc
end

% Perform reconstruction
cd (path_sd3d);
files_jac = dir('*.jac');
files_invjac = dir('*.invjac');
files_rmap = dir('*.rmap');
cd(path_prepro);
files_prepro = dir('*gsr.prepro');

for nsub = 1:length(files_jac)
    
    tic;
    rmapFileName = [path_sd3d filesep files_rmap(nsub).name];
    invjacFileName = [path_sd3d filesep files_invjac(nsub).name];
    preproFileName = [path_prepro filesep files_prepro(nsub).name];
    %  Don't save volume files as they are not needed and are very big
    [dotimg, dotimgFileName] = DOTHUB_reconstruction(preproFileName,[],invjacFileName,rmapFileName,'saveVolumeImages',false);
    toc
    
end





