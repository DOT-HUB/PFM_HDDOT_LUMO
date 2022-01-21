
function [Yparcels_HbO, Yparcels_HbR]= ParcellateNodes(dotimg, gmSurfaceMeshLabels)

% Function to parcellate a grey matter surface mesh into anatomical parcels. 

% Input: 
% dotimg                dotimg file as per DOT-HUB reconstruction pipeline
% gmSurfaceMeshLabels   gmSurfaceMesh file that has labels defining which
%                       node corresponds to which anatomical region

% Output:               
% Yparcels_HbO,HbR      matrices of size nFramesxnParcels indicating which
%                       node is in which anatomical parcels. nParcels will
%                       be the same as the number of labels
  

[a,~,c] = unique(gmSurfaceMeshLabels(:,1));
nParcels = length(a);
nNodes = length(c);
a(1) % Check if it is zeros
Ytrans = zeros(nParcels, nNodes);

for i = 1:nParcels
    
    tmpLog = (c==i);
    Ytrans(i,tmpLog) = 1/sum(tmpLog);
    
end

Y_parcelstmpHbO = Ytrans*dotimg.hbo.gm';
Yparcels_HbO = Y_parcelstmpHbO';

Y_parcelstmpHbR = Ytrans*dotimg.hbr.gm';
Yparcels_HbR = Y_parcelstmpHbR';


% % ABOVE IT IS COMPUTING THE MEAN, EQUIVALENT TO BELOW
% yparcel_hbo = zeros(size(dotimg.hbo.gm,1),nParcels);
% yparcel_hbr = zeros(size(dotimg.hbr.gm,1),nParcels);
% 
% for i = 1:nParcels
%     tmpLog = (c==i);   
%     yparcel_hbo(:,i) = mean(dotimg.hbo.gm(:, tmpLog==1),2);
%     yparcel_hbr(:,i) = mean(dotimg.hbr.gm(:, tmpLog==1),2);    
% end
    

end

