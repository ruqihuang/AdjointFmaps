clc; clear; close all; 
addpath(genpath('data')); 
addpath(genpath('external')); 
addpath(genpath('utils')); 

load human_poses.mat

nshapes = length(shapes); 
kEig = 60; 
L = cell(nshapes, 1); 
for i = 1:nshapes
    L{i} = compute_laplacian_basis(shapes{i}, kEig); 
end

% compute ground-truth funcitonal maps
C = cell(nshapes);
for i = 1:nshapes
    for j = 1:nshapes
        if i ~= j
            C{i, j} = L{j}.evecs'*L{j}.A*L{i}.evecs; 
        else
            C{i, j} = eye(kEig); 
        end
    end
end

U= extract_latent_basis(C); 

l = 40; 
B = U.evecs(:, 1:l); 

X = cell(size(C)); 
for i = 1:length(C)
    for j = 1:length(C)
        X{i, j} = C{j, i}';
    end
end

W = extract_latent_basis(X); 


E = B'*W.mat*B; 
E = (E+E')/2; 

[u, v] = eig(E); 
[v, idx] = sort(diag(v), 'descend'); 
u = u(:, idx); 

figure; 
for i = 1:5
    for j = 1:5
        if j == 1
            subplot(5, 5, (i-1)*5+j); plotMesh(L{j}, (L{j}.evecs*U.bases{j}(:, 1:l)*u(:, i)).^2, -90, -6); 
        elseif j == 4
            subplot(5, 5, (i-1)*5+j); plotMesh(L{j}, (L{j}.evecs*U.bases{j}(:, 1:l)*u(:, i)).^2, -180, -6); 
        else
            subplot(5, 5, (i-1)*5+j); plotMesh(L{j}, (L{j}.evecs*U.bases{j}(:, 1:l)*u(:, i)).^2, 90, 6); 
        end
    end
end
% print('F7L', '-dpng', '-r400'); 
        
        