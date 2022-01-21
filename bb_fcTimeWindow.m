function [group_hbo, group_hbr, tw_order] = bb_fcTimeWindow(data_hbo, data_hbr, tw_vect, mask, groundT_hbo, groundT_hbr)

[nrows, ncols] = size(data_hbo);
windowstep = tw_vect(1)/2; % half tw duration, step between each window. If < windowlength, windows overlap

idx_tw = 1;
group_hbo = [];
group_hbr = [];

for i = 1:length(tw_vect)
    %inputs:
    windowlength = tw_vect(i);  %size of windows
    
    % Create windows
    %switch signals along 3rd dimension:
    seg_data_hbo = data_hbo(bsxfun(@plus, ...
        bsxfun(@plus, ...
        1:windowlength, ...
        (0:windowstep:ncols-windowlength).') * nrows, ...
        permute(1-nrows:0, [1 3 2])));
    
    seg_data_hbr = data_hbr(bsxfun(@plus, ...
        bsxfun(@plus, ...
        1:windowlength, ...
        (0:windowstep:ncols-windowlength).') * nrows, ...
        permute(1-nrows:0, [1 3 2])));
    
    % Compute graph metrics
    % For each time window
    % Store values
    for j = 1:size(seg_data_hbo,1)
        
        % Fisher transformed for normality
        temp_hbo = atanh(corr(squeeze(seg_data_hbo(j,:,:))));
        temp_hbo = temp_hbo(mask);
        rsquared_hbo = corr(temp_hbo, groundT_hbo)^2;
        
        group_hbo(1, idx_tw) = i;
        group_hbo(2, idx_tw) = rsquared_hbo;
        
        % Fisher transformed for normality
        temp_hbr = atanh(corr(squeeze(seg_data_hbr(j,:,:))));
        temp_hbr = temp_hbr(mask);
        rsquared_hbr = corr(temp_hbr, groundT_hbr)^2;
        
        group_hbr(1, idx_tw) = i;
        group_hbr(2, idx_tw) = rsquared_hbr;
        
        idx_tw = idx_tw + 1;
        
    end
   clc
 
end

tw_order = group_hbo(1,:);
group_hbo = group_hbo(2,:);                     
group_hbr = group_hbr(2,:);
                  
                         
                       

