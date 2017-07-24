%In the end, we reproduce figure(1) as Figure 6 and figure(2)
%as Figure 13 in the paper. 

clc; clear; close all; 
addpath(genpath('data')); 
addpath(genpath('external')); 
addpath(genpath('utils')); 

S1 = read_off_shape('1.off'); 
S2 = read_off_shape('11.off'); 

kEig = 60; 
L{1} = compute_laplacian_basis(S1, kEig); 
L{2} = compute_laplacian_basis(S2, kEig); 

% compute ground-truth funcitonal maps
C = cell(2);
for i = 1:2
    for j = 1:2
        if i ~= j
            C{i, j} = L{j}.evecs'*L{j}.A*L{i}.evecs; 
        else
            C{i, j} = eye(kEig); 
        end
    end
end

% compute the shape difference operators
[u1, v1] = eig(C{1, 2}'*C{1, 2}); 
[v1, ind] = sort(diag(v1), 'descend'); 
u1 = u1(:, ind); 
f = zeros(12500, 1); 
for i = 1:5
    f = f + v1(i)*(L{1}.evecs*u1(:, i)).^2; 
end



[u2, v2] = eig(C{2, 1}'*C{2, 1});
[v2, ind] = sort(diag(v2), 'descend'); 
u2 = u2(:, ind); 
g = zeros(12500, 1); 
for i = 1:5
    g = g+v2(i)*(L{2}.evecs*u2(:, i)).^2; 
end


U = extract_latent_basis(C); 
l = 49; 

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
h1 = zeros(12500, 1); 
for i = 1:5
    h1 = h1 + v(i)*(L{1}.evecs*(U.bases{1}(:, 1:l)*u(:, i))).^2; 
end

h2 = zeros(12500, 1); 
for i = 1:5
    h2 = h2 + v(i)*(L{2}.evecs*(U.bases{2}(:, 1:l)*u(:, i))).^2; 
end


figure(1); 
subplot(2, 4, 1); 
plotMesh(L{1}, f, -90, -6);
subplot(2, 4, 5); 
plotMesh(L{1}, f, 90, 6); 
subplot(2, 4, 2); 
plotMesh(L{2}, g, 90, 6); 
subplot(2, 4, 6); 
plotMesh(L{2}, g, -90, -6); 
subplot(2, 4, 3); 
plotMesh(L{1}, h1, -90, -6); 
subplot(2, 4, 4); 
plotMesh(L{2}, h2, 90, 6); 
subplot(2, 4, 7); 
plotMesh(L{1}, h1, 90, 6); 
subplot(2, 4, 8); 
plotMesh(L{2}, h2, -90, -6); 

figure(2); 
for i = 1:5
    subplot(2, 10, i); plotMesh(L{1}, (L{1}.evecs*u1(:, i)).^2, -90, -6);
    subplot(2, 10, i+5); plotMesh(L{2}, (L{2}.evecs*u2(:, i)).^2, 90, 6); 
end
subplot(2, 10, 2); plotMesh(L{1}, (L{1}.evecs*u1(:, 2)).^2, 90, 6);
subplot(2, 10, 6); plotMesh(L{2}, (L{2}.evecs*u2(:, 1)).^2, -90, -6);
subplot(2, 10, 8); plotMesh(L{2}, (L{2}.evecs*u2(:, 3)).^2, -90, -6);

for i = 1:10
    subplot(2, 10, 10+i); plotMesh(L{1}, (L{1}.evecs*(U.bases{1}(:, 1:l)*u(:, i))).^2, -90, -6); 
end
for i = [3, 4, 6]
    subplot(2, 10, 10+i); plotMesh(L{1}, (L{1}.evecs*(U.bases{1}(:, 1:l)*u(:, i))).^2, 90, 6);
end



