
clear; close all; clc

path_dotimg = ('/Volumes/CAM_data/neoLAB/dotimg/filter');
cd(path_dotimg);
files_prep = dir('*.mat');
groupFC_z = [];

for nsub = 1:length(files_prep)
    
    cd(path_dotimg);
    load(files_prep(nsub).name, '-mat');
    groupFC_z = cat(3, groupFC_z, sub_data.FC_zparc);
    
end

fig1 = figure; set(fig1,'units', 'normalized', 'outerposition', [0 0 0.8 0.8], 'Color', [1 1 1]);
FC_group = tanh(mean(groupFC_z, 3));
imagesc(FC_group, [-1 1]); colormap jet
set(gca, 'FontSize', 24)
xlabel('Parcels  HbO-HbR')
ylabel('Parcels  HbO-HbR')