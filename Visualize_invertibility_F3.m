clear all; close all; clc; 

addpath(genpath('data')); 
addpath(genpath('external')); 
addpath(genpath('utils')); 

maps = load('Faust_030_090_maps');

%%
srcmesh = 'tr_reg_030.off';
tarmesh = 'tr_reg_090.off';

Src = read_off_shape(['data/FAUST/' srcmesh]);
fprintf('%s vs %s\n',srcmesh,tarmesh);
%fprintf('reading the source shape...');tic;
Src = compute_laplacian_basis(Src, 1);
%fprintf('done (found %d vertices)\n',Src.nv);toc;

%fprintf('reading the target shape...');tic;
Tar = read_off_shape(['data/FAUST/' tarmesh]);
Tar = compute_laplacian_basis(Tar, 1);

%%
map_Ad = maps.map_Ad;
load('MyColormaps');
FigHandle = figure('Position', [100, 100, 800, 450]);
colormap(mycmap);
 f1 = ones(size(map_Ad));
 f1(1000) = 0;
 plot_function_faust(Src,f1);
 xlabelm('Source Shape',-0.3);
%  saveas(FigHandle,'source_shape.png');

% title('Source Shape','FontSize',16);
hold off;
 %%
FigHandle = figure('Position', [100, 100, 800, 450]);
hold on;
 colormap(mycmap);

Tar.Av = full(diag(Tar.A));
Tar.As = sum(Tar.Av);
subplot(1,3,1);
[m1,m2] = ismapped(maps.map_NA);
plot_function_faust(Tar,m1);
title(sprintf('Coverage = %.1f%%\n',100*sum(Tar.Av(m2))/Tar.As),'FontSize',16);
xlabelm('Regular Fmap',-0.3);

subplot(1,3,2);
[m1,m2] = ismapped(maps.map_Cf);
plot_function_faust(Tar,m1);
title(sprintf('Coverage = %.1f%%\n',100*sum(Tar.Av(m2))/Tar.As),'FontSize',16);
xlabelm('[ERGB]',-0.15);

subplot(1,3,3);
[m1,m2] = ismapped(maps.map_Ad);
plot_function_faust(Tar,m1);
title(sprintf('Coverage = %.1f%%\n',100*sum(Tar.Av(m2))/Tar.As),'FontSize',16);
xlabelm('Adjoint Regularization',-0.40);
% saveas(FigHandle,'coverage_noICP.png');

%%
FigHandle = figure('Position', [100, 100, 800, 450]);
hold on;
 colormap(mycmap);

subplot(1,3,1);
[m1,m2] = ismapped(maps.map_NAICP);
plot_function_faust(Tar,m1);
xlabelm('Regular Fmap + ICP',-0.5);
title(sprintf('Coverage = %.1f%%\n',100*sum(Tar.Av(m2))/Tar.As),'FontSize',20);


subplot(1,3,2);
[m1,m2] = ismapped(maps.map_CfICP);
plot_function_faust(Tar,m1);
title(sprintf('Coverage = %.1f%%\n',100*sum(Tar.Av(m2))/Tar.As),'FontSize',20);
xlabelm('[ERGB] + ICP',-0.25);


subplot(1,3,3);
[m1,m2] = ismapped(maps.map_AdICP);
plot_function_faust(Tar,m1);
title(sprintf('Coverage = %.1f%%\n',100*sum(Tar.Av(m2))/Tar.As),'FontSize',20);
xlabelm('Adjoint Regularization + ICP',-0.5);
% saveas(FigHandle,'coverage_wicp.png');



%%


%%
% samples = dijkstra_fps(Src,200);
% SS = Src;
% TS = Tar;
% SS.X = (roty(30)*SS.Pts')';
% TS.X = (roty(30)*TS.Pts')';
%TS.X(:,3) = SS.X(:,1)+10;

% figure(4);
% visualize_map_lines(SS,TS,maps.map_NAICP,samples);
% title('[Regular Fmap] + ICP','FontSize',20,'FontWeight','b');
% axis off;
% view([-2 90]);
% 
% figure(5);
% visualize_map_lines(SS,TS,maps.map_AdICP,samples);
% title('[Adjoint Regularization] + ICP','FontSize',20,'FontWeight','b');
% axis off;
% view([-2 90]);


%saveas(gcf,'ergb_mapvis.png');
%visualize_map_lines(Src,Tar,map_NAICP);

