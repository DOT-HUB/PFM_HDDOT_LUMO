clear; close all; clc
path_data = '/Volumes/CAM_data/neoLAB/dotimg/ssr';
cd(path_data)
sub = dir('*ssr.mat');
path_sub ='/Volumes/CAM_data/neoLAB/group/Data4DualRegression';


for nsub = 1:length(sub)
    clc
    nsub
    cd(path_data)
    load(sub(nsub).name, '-mat');
    hbo = sub_data.hbo_zparc;
    hbr = sub_data.hbr_zparc;
    
    cd(path_sub)
    save([sub(nsub).name(1:end-7) 'hbo.mat'],'hbo')
    save([sub(nsub).name(1:end-7) 'hbr.mat'],'hbr')

end
    
    