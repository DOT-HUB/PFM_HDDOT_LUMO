
% CONSIDER null_model_und_sign for creating
% networks to test null hypothesis
clear all; close all; clc
addpath(genpath('/Users/borjablanco/Documents/GitHub/Hyper24RSFC_scripts/scripts/BB_scripts'))
load('/Volumes/CAM_data/neoLAB/group/test_graphT')

% ------- Modularity
[M,Q] = modularity_und((test>0.2).*test,1);
[ML,QL] = community_louvain((test>0).*test,1, [], 'modularity')

% ------- Global efficiency
% Only works with positive values
E = efficiency_wei((test>0.6).*test,0);

% ------- Clustering coefficient
% Only works with positive values
C = mean(clustering_coef_wu((test>0).*test));

% ------- Assortativity
A = assortativity_wei((test>0).*test, 0)

% Transitivity
%T = transitivity_wu();

% ------- Participation coefficient
% Needs a vector assigning community to each parcel eg.g Ci from modularity
Pneg = participation_coef(test, parc_net(:,2), 0);

% ------- Betweenes centrality
% needs distance matrix
test_inv = find(test); test(test_inv)=1./test(test_inv);% invert weights
N = 154; % number of nodes
BC = betweenness_wei(test_inv);

% normalize between 0 and 1 by the number of parcels
BCn = BC/((N-1)*(N-2));

% ------- Rich club coefficient
RC = rich_club_wu ((test>0).*test);

% ------- Degree
% Apply threshold mask (e.g., 10% higher), otherwise all ones
D = degrees_und(test);







