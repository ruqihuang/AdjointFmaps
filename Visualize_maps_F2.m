clc; clear; close all; 

addpath(genpath('data')); 
addpath(genpath('external')); 
addpath(genpath('utils')); 

maps = load('Faust_030_090_maps.mat');

%%
srcmesh = 'tr_reg_030.off';
tarmesh = 'tr_reg_090.off';

Src = read_off_shape(['data/FAUST/' srcmesh]);
fprintf('%s vs %s\n',srcmesh,tarmesh);
fprintf('reading the shapes...');tic;
Src = compute_laplacian_basis(Src, 1);
Tar = read_off_shape(['data/FAUST/' tarmesh]);
Tar = compute_laplacian_basis(Tar, 1);
fprintf('done\n');

%%
load('MyColormaps');

FigHandle1 = figure('Position', [100, 100, 800, 600]);
colormap(mycmap);

%%
fprintf('sampling the source points...');
samples = dijkstra_fps(Src,200);
fprintf('done\n');

SS = Src;
TS = Tar;
SS.X = (roty(30)*SS.Pts')';
TS.X = (roty(30)*TS.Pts')';

%%
visualize_map_lines(SS,TS,maps.map_NAICP,samples);
title('[Fmap] + ICP','FontSize',20,'FontWeight','b');
axis off;
view([-2 90]);

FigHandle2 = figure('Position', [200, 200, 800, 600]);
visualize_map_lines(SS,TS,maps.map_AdICP,samples);
title('[Adjoint] + ICP','FontSize',20,'FontWeight','b');
axis off;
view([-2 90]);
