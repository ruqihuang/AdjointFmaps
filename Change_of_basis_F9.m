close all; clear; clc;
addpath(genpath('data')); 
addpath(genpath('external')); 
addpath(genpath('utils')); 

%% load shapes, ground-truth p2p map and the recovery result from [RMC15]. 

% Choose different pair of shapes by un-commenting one of the following: 

% elephant vs. horse
load data_elephant_horse_subsampled_matches.mat;
load elephanthorse_RMC15.mat
S1.surface = shapes{1};
S2.surface = shapes{2}; 
gt_match = 1:length(S1.surface.X);
gt_match1 = gt_match; 
nk = 60; 
nk2 = 40;

% % conformal bunny
% load bunny_RMC15.mat
% S1 = read_off_shape('flow00.off');
% S2 = read_off_shape('flow50.off');
% gt_match = 1:length(S1.surface.X); 
% gt_match1 = gt_match; 
% nk = 60; 
% nk2 = 30; 

% 
% % cat vs. lion
% load catlion_RMC15.mat
% S1 = read_off_shape('cat-00.off');
% S2 = read_off_shape('lion-00.off');
% load cat-0_lion-0_rough_map.mat
% gt_match = nnidx; 
% gt_match1 = innidx; 
% nk = 60;
% nk2 = 30;

% % deformed faces
% load faces_RMC15.mat
% S1 = read_off_shape('hao_li_001.off'); 
% S2 = read_off_shape('deformed_hao_li_003.off'); 
% gt_match = 1:length(S1.surface.X); 
% gt_match1 = gt_match; 
% nk = 60;
% nk2 = 30;
% 


% Compute the LB eigenbasis
S1 = compute_laplacian_basis(S1, nk);
S2 = compute_laplacian_basis(S2, nk);


C12 = S2.evecs'*S2.A*S1.evecs(gt_match1, :); 
e = ones(length(S2.evecs), 1); 
e = S2.evecs'*S2.A*e; 
k1_const = S1.evecs*C12'*e;
K = S1.evecs'*S1.A*diag(k1_const)*S1.evecs;
X12 = C12/K;

cvx_begin quiet
    variable H(nk, nk) semidefinite
    minimize( norm(X12 - H*C12, 'fro'))        
cvx_end

[e,v] = eig(diag(S2.evals),H);
e = e*diag(1./sqrt(diag(e'*H*e)));
[v, order] = sort(diag(v));
e = e(:,order);

%%
% Compute the functional map with respect to this new basis.
C2 = (S2.evecs*e)\S1.evecs(gt_match1, 1:nk);

C2r = C2(1:nk,1:nk2);


%% Convert p2p maps
[~, matches_nn] = icp_refine(S1.evecs, S2.evecs, C12, 0); 
[~, matches_cb] = icp_refine(S1.evecs,(S2.evecs*e),C2r,0);

error_nn = eval_geodesic_error(S1, matches_nn, gt_match1)/sqrt(6); 
error_rmc = eval_geodesic_error(S1, matches3, gt_match1)/sqrt(6); 
error_cb = eval_geodesic_error(S1, matches_cb, gt_match1)/sqrt(6); 


figure(1);    
nr = length(reshape(error_nn,[],1));
plot(sort(reshape(error_nn,[],1))*100,linspace(0,1,nr)*100,'-k','LineWidth',2);
hold on; 
plot(sort(reshape(error_rmc,[],1))*100,linspace(0,1,nr)*100,'-b','LineWidth',2);
hold on;
plot(sort(reshape(error_cb,[],1))*100,linspace(0,1,nr)*100,'-r','LineWidth',2);
xlim([0, 25]);
ylim([0, 100]); 
ylabel('% of Correspondences', 'FontSize',25);
xlabel('% of Geodesic Error', 'FontSize',25);
leg = legend('NN-search', '[RMC15]', 'Our scheme' ,'Location', 'southeast'); 
set(leg, 'FontSize', 25); 
title('Elephant vs horse', 'FontSize', 20); 

figure(2); 
subplot(1, 3, 1); plotMesh(S1, zeros(length(S1.surface.X), 1), 0, 90, 'f'); 
subplot(1, 3, 2); plotMesh(S2, error_rmc, 0, 90, 'f'); 
subplot(1, 3, 3); plotMesh(S2, error_cb, 0, 90, 'f'); 



