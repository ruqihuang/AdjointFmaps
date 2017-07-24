clc; clear; close all; 
addpath(genpath('data')); 
addpath(genpath('external')); 
addpath(genpath('utils')); 

load galloping_horses.mat

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

[u, v] = eig(E); 
[v, idx] = sort(diag(v), 'descend'); 
u = u(:, idx); 

figure; 
for i = 1:5
    for j = 1:nshapes
        if (i == 2) || (i == 5)
            subplot(5, nshapes, (i-1)*nshapes+j); plotMesh(L{5-j}, (L{5-j}.evecs*(U.bases{5-j}(:, 1:l)*u(:, i))).^2, 90, -90); 
        else
            subplot(5, nshapes, (i-1)*nshapes+j); plotMesh(L{5-j}, (L{5-j}.evecs*(U.bases{5-j}(:, 1:l)*u(:, i))).^2, -90, 90); 
        end
    end
end
% print('F7R', '-dpng', '-r400'); 
        
        