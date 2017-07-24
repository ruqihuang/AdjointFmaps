% Readme: we first select randomly 25 pairs of deformed spheres, and then
% recover the p2p maps from the ground-truth functional maps. A fixed set
% of random pairs and pre-computed result from the paper [RMC15] are loaded
% in line 19~20. (Unfortunately we lost the random pairs
% used in producing the result in the paper) 

clear; close all; clc; 
addpath(genpath('data')); 
addpath(genpath('external')); 
addpath(genpath('utils')); 

% %% generate random pairs of deformed spheres. 
% ra = randi(64, 50, 2); 
% ra = unique(ra, 'rows'); 
% idx = find(ra(:, 1) - ra(:, 2) ~= 0); 
% ra = ra(idx, :); 
% t = randperm(33); 
% ra = ra(t(1:25), :); 

%% Pre-load random pairs and results from [RMC15]
load rand_sph.mat;
load cpd_results.mat; 

%% Recovery
E1 = zeros(1922, 25); 
E2 = E1; 
E3 = E1; 
for i = 1:25
    test_shape1 = sprintf('deform_sph_%d.off', ra(i, 1)); 
    test_shape2 = sprintf('deform_sph_%d.off', ra(i, 2)); 

    S1 = read_off_shape(test_shape1);
    S2 = read_off_shape(test_shape2);


    % Number of basis functions
    nk = 60;
    % dimensionality of the functional map (nk x nk2)
    nk2 = 20;

    % Compute the LB eigenbasis
    S1 = compute_laplacian_basis(S1, nk);
    S2 = compute_laplacian_basis(S2, nk);


    % Compute initial (slightly perturbed identity) map
    P2 = 1:length(S1.evecs); 
    P1 = 1:length(S2.evecs); 
    C = S2.evecs'*S2.A*S1.evecs(P1, :); 


    C12 = C;
    e = ones(length(S2.evecs), 1); e = S2.evecs'*S2.A*e; 
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

    % Compute the functional map with respect to this new basis.
    C2 = (S2.evecs*e)\S1.evecs(P1, 1:nk);
    C2r = C2(1:nk,1:nk2);


    % Covert to p2p map
    [~, matches1] = icp_refine(S1.evecs, S2.evecs, C, 0); 
    [~, matches2] = icp_refine(S1.evecs,(S2.evecs*e),C2r,0);

    error1 = eval_geodesic_error(S1, matches1, P1)/sqrt(6); 
    error2 = eval_geodesic_error(S1, M(:, i), P1)/sqrt(6); 
    error3 = eval_geodesic_error(S1, matches2, P1)/sqrt(6);

    E1(:, i) = error1; 
    E2(:, i) = error2; 
    E3(:, i) = error3; 

    fprintf(sprintf('the %d-th pair done.\n', i)); 
end

%% Plot the results. 
nr = length(reshape(E1, [], 1));
figure; plot(sort(reshape(E1, [], 1))*100, linspace(0, 1, nr)*100, 'k', 'LineWidth', 3); 
hold on; plot(sort(reshape(E2, [], 1))*100, linspace(0, 1, nr)*100, 'b', 'LineWidth', 3); 
hold on; plot(sort(reshape(E3, [], 1))*100, linspace(0, 1, nr)*100, 'r', 'LineWidth', 3); 
ylabel('% of Correspondences', 'FontSize',25);
xlabel('% of Geodesic Error', 'FontSize',25);
ylim([80, 100]); 
xlim([0, 8]); 
leg = legend('NN-search', '[RMC15]', 'Our scheme' ,'Location', 'southeast'); 
set(leg, 'FontSize', 25);
% print('sphere_collection', '-r400', '-dpng'); 
