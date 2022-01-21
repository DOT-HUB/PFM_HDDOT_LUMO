% parrcel sensitivity


% Function find the parrcel sensitivity for the group level GMSensitivity Mask 

% Input: 
% GMSensitivty          binarized GM sensitivity mask at the group leve
% gmSurfaceMeshLabels   gmSurfaceMesh file that has labels defining which
%                       node corresponds to which anatomical region


[a,~,c]=unique(gmSurfaceMeshLabels(:,1));
nParcels=length(a);
nNodes=length(c);

Ytrans=zeros(nParcels, nNodes);


for i=1:nParcels
    tmpLog=c==i;
    Ytrans(i,tmpLog)=1/sum(tmpLog);
end


Y_parcelstmpMask=Ytrans*GMmaskSumBinary'
Yparcels_mask=Y_parcelstmpMask';


