clc; clear; close all; 
addpath(genpath('data')); 
addpath(genpath('external')); 
addpath(genpath('utils')); 

S1 = read_off_shape('deform_sph_38.off'); 
S2 = read_off_shape('deform_sph_8.off'); 

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

[u2, v2] = eig(C{2, 1}'*C{2, 1});
[v2, ind] = sort(diag(v2), 'descend'); 
u2 = u2(:, ind); 

s = max(max(v1), max(v2)); 


U= extract_latent_basis(C); 



% Select the consistent basis
l = 52; 
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


f = u(:, 1); 
f1 = L{1}.evecs*(U.bases{1}(:, 1:l)*f); 
f2 = L{2}.evecs*(U.bases{2}(:, 1:l)*f); 
g = u(:, 4); 
g1 = L{1}.evecs*(U.bases{1}(:, 1:l)*g); 
g2 = L{2}.evecs*(U.bases{2}(:, 1:l)*g); 

figure; 
subplot(3, 4, 1); 
plot(v1, 'LineWidth', 2); 
xlabel('Index', 'FontWeight','bold','FontSize',30);
ylabel('Spectrum of D_{M, N}', 'FontWeight','bold', 'FontSize',30);
ylim([0, s+0.1]); 
xlim([0, 60]); 
set(gca,'fontsize',5)
hold on; 
plot(1,v1(1),'r.', 'MarkerSize',10); 
hold on; 
plot(60, v1(60), 'k.', 'MarkerSize',10); 

subplot(3, 4, 2); 
plot(v2, 'LineWidth', 2); 
xlabel('Index', 'FontWeight','bold','FontSize',30);
ylabel('Spectrum of D_{N, M}', 'FontWeight','bold', 'FontSize',30);
ylim([0, s+0.1]); 
xlim([0, 60]); 
set(gca,'fontsize',5)
hold on; 
plot(1,v2(1),'r.', 'MarkerSize',10); 
hold on; 
plot(60, v2(60), 'k.', 'MarkerSize',10); 

subplot(3, 4, 3); 
plot(v, 'LineWidth', 2); 
xlabel('Index', 'FontWeight','bold','FontSize',30);
ylabel('Spectrum of V^X', 'FontWeight','bold', 'FontSize',30);
xlim([0, 60])
set(gca,'fontsize',5)
hold on; 
plot(1,v(1), 'r.', 'MarkerSize',10); 
hold on; 
plot(4, v(4), 'k.', 'MarkerSize',10); 

subplot(3, 4, 5); 
plotMesh(L{1}, (L{1}.evecs*u1(:, 1)).^2, 180, -60);

subplot(3, 4, 6); 
plotMesh(L{2}, (L{2}.evecs*u2(:, 1)).^2, 180, -60);

subplot(3, 4, 7); 
plotMesh(L{1}, f1.^2, 180, -60); 

subplot(3, 4, 8); 
plotMesh(L{2}, f2.^2, 180, -60); 

subplot(3, 4, 9); 
plotMesh(L{1}, (L{1}.evecs*u1(:, 60)).^2, 180, -60); 

subplot(3, 4, 10); 
plotMesh(L{2}, (L{2}.evecs*u2(:, 60)).^2, 180, -60);

subplot(3, 4, 11); 
plotMesh(L{1}, g1.^2, 180, -60); 

subplot(3, 4, 12); 
plotMesh(L{2}, g2.^2, 180, -60); 
