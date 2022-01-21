function [gmetrics_hbo, gmetrics_hbr, tw_order] = bb_graphTimeWindow(data_hbo, data_hbr, tw_vect)

% Compute data size
[nrows, ncols] = size(data_hbo);

% half tw(1) duration = 30s
% This is the step between each window.
windowstep = tw_vect(1)/2; 

idx_tw = 1;
gmetrics_hbo = [];
gmetrics_hbr = [];
tw_order = [];

% Compute metrics for a set of density thresholds
dens = 1:1:50;

for i = 1:length(tw_vect)
    %inputs:
    windowlength = tw_vect(i);  % window size (1, 2, 3 ... minutes)
    
    % Create segmented data for each time window
    % switch signals along 3rd dimension:
    % The number of time windows depends on windowlength
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
    % For each time window and for each density threshold
    % Store values
    for j = 1:size(seg_data_hbo,1)
        
        % Select data for the time window
        temp_hbo = corr(squeeze(seg_data_hbo(j,:,:)));
        temp_hbr = corr(squeeze(seg_data_hbr(j,:,:)));
        
        % Compute graph metrics at different density thresholds
        for th = 1:length(dens)
            
            % Calculate threshold for density values, e.g. 5%
            vect_hbo = triu(temp_hbo,1); % vectorize FC matrix
            vect_hbo = sort(vect_hbo(vect_hbo~=0), 'descend');
            th_hbo = vect_hbo(round(dens(th)/100*length(vect_hbo))+1);
            
            vect_hbr = triu(temp_hbr,1);
            vect_hbr = sort(vect_hbr(vect_hbr~=0), 'descend');
            th_hbr = vect_hbr(round(dens(th)/100*length(vect_hbr))+1);
            
            % Apply proportional threshold to fc matrix
            ses_hbo = (temp_hbo>=th_hbo).*temp_hbo;
            ses_hbr = (temp_hbr>=th_hbr).*temp_hbr;
            
            % Remove NaN from main diag
            ses_hbo(isnan(ses_hbo)) = 1;
            ses_hbr(isnan(ses_hbr)) = 1;
            
            % Compute graph metrics
            % Global efficiency
            ge_hbo = efficiency_wei((ses_hbo>0).*ses_hbo,0);
            ge_hbr = efficiency_wei((ses_hbr>0).*ses_hbr,0);
            
            % Modularity
            [~, cl_hbo] = community_louvain(ses_hbo,1, [], 'negative_asym');
            [~, cl_hbr] = community_louvain(ses_hbr,1, [], 'negative_asym');
            
            [~, ml_hbo] = community_louvain((ses_hbo>0).*ses_hbo,1, [], 'modularity');
            [~, ml_hbr] = community_louvain((ses_hbr>0).*ses_hbr,1, [], 'modularity');
            
            % Clustering coefficient
            cc_hbo = mean(clustering_coef_wu((ses_hbo>0).*ses_hbo));
            cc_hbr = mean(clustering_coef_wu((ses_hbr>0).*ses_hbr));

            % Assortativity
            as_hbo = assortativity_wei((ses_hbo>0).*ses_hbo, 0);
            as_hbr = assortativity_wei((ses_hbr>0).*ses_hbr, 0);

            % Store results for each tw and density threshold
            gmetrics_hbo(1, idx_tw, th) = ge_hbo;
            gmetrics_hbo(2, idx_tw, th) = cl_hbo;
            gmetrics_hbo(3, idx_tw, th) = ml_hbo;
            gmetrics_hbo(4, idx_tw, th) = cc_hbo;
            gmetrics_hbo(5, idx_tw, th) = as_hbo;
           
            gmetrics_hbr(1, idx_tw, th) = ge_hbr;
            gmetrics_hbr(2, idx_tw, th) = cl_hbr;
            gmetrics_hbr(3, idx_tw, th) = ml_hbr;
            gmetrics_hbr(4, idx_tw, th) = cc_hbr;
            gmetrics_hbr(5, idx_tw, th) = as_hbr;
            
        end
        tw_order(idx_tw) = i;      
        idx_tw = idx_tw + 1;
    end
        
end









