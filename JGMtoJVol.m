function [J_vol_wav1,J_vol_wav2, J_GM_wav1, J_GM_wav2] = JGMtoJVol(jacFileName,rmapFileName,origMeshFileName)

%written by Liam Collins-Jones, UCL, 2020
% updated by Julie Uchitel, Jan 2021
load(origMeshFileName,'-mat','gmSurfaceMesh');
gmSurfaceMeshCommon = gmSurfaceMesh;
load(rmapFileName,'-mat','headVolumeMesh','gmSurfaceMesh','vol2gm');
basis = [30 30 30];
fineBasis = basis.*2;
disp('Mapping basis Jacobian to volume');
eltp = ones(length(headVolumeMesh.elem),1)*3;
hMesh = toastMesh(headVolumeMesh.node(:,1:3),headVolumeMesh.elem(:,1:4),eltp);
hBasis = toastBasis(hMesh,basis,fineBasis);
load(jacFileName,'-mat','J');
J_GM_wav1 = J{1,1}.gm;
J_GM_wav2 = J{2,1}.gm;
J_basis_wav1 = J{1,1}.basis;
J_basis_wav2 = J{2,1}.basis;
disp('Mapping basis Jacobian to volume');
%Map back to volume
J_vol_wav1 = zeros(size(J_basis_wav1,1),length(headVolumeMesh.node));
J_vol_wav2 = zeros(size(J_basis_wav1,1),length(headVolumeMesh.node));
for i = 1:size(J_basis_wav1,1);
    i
    J_vol_wav1(i,:) = hBasis.Map('S->M',J_basis_wav1(i,:));
    J_vol_wav2(i,:) = hBasis.Map('S->M',J_basis_wav2(i,:));
end
