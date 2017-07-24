% This code implements a basic version of the algorithm described in:
%
% Informative Descriptor Preservation via Commutativity for Shape Matching,
% Dorian Nogneng and Maks Ovsjanikov, Proc. Eurographics 2017
%
% To try it, simply run this file in MATLAB. This should produce
% a map (correspondence) between a pair of meshes from the FAUST dataset,
% and create an image that visualizes this correspondence.
%
% This code was written by Etienne Corman and modified by Maks Ovsjanikov.

%clear all; close all;

clear; close all;
addpath(genpath('data')); 
addpath(genpath('external')); 
addpath(genpath('utils')); 
load Tosca_ERGB.mat; 

%% Load meshes and compute Laplacian eigendecomposition

% Number of basis vectors for computing the functional map.
% Larger is usually better (more accurate results) but somewhat slower.
numEigsSrc = 60;
numEigsTar = 60;

meshes = dir('data/TOSCA_simplified/*.off');

rng('default');
rng(2);

categories = {'cat','centaur','david','dog','horse','michael','victoria','wolf'};

tars = [];
for i=1:length(categories)
    a = dir(['data/TOSCA_simplified/' categories{i} '*.off']);
    tars = [tars; randperm(length(a))'+length(tars)];
    fprintf('%s %d\n',categories{i},length(a));
end

epss1 = 500;
epss = 0.000000453999298;

allerrsA = [];
allerrsICPA = [];

allerrsNA = [];
allerrsICPNA = [];

meanerrs = [];
for k=1:1:length(meshes)
    srcmesh = meshes(k).name;
    tarmesh = meshes(tars(k)).name;
    
    Src = read_off_shape(['data/TOSCA_simplified/' srcmesh]);
    fprintf('%s vs %s\n',srcmesh,tarmesh);
    %fprintf('reading the source shape...');tic;
    Src = compute_laplacian_basis(Src, 200);
    %fprintf('done (found %d vertices)\n',Src.nv);toc;
    
    %fprintf('reading the target shape...');tic;
    Tar = read_off_shape(['data/TOSCA_simplified/' tarmesh]);
    Tar = compute_laplacian_basis(Tar, 200);
    %fprintf('done (found %d vertices)\n',Tar.nv);toc;
    
    % a few landmark correspondences (to avoid symmetry flipping).
    landmarks1 = (500:1000:3000)';
    %landmarks1 = dijkstra_fps(Src,20);
    landmarks2 = landmarks1;
    
    landmarks = [landmarks1 landmarks2(:,1)];
    
    SrcLaplaceBasis = Src.evecs; SrcEigenvalues = Src.evals;
    TarLaplaceBasis = Tar.evecs; TarEigenvalues = Tar.evals;
    Src.evecs = SrcLaplaceBasis(:,1:numEigsSrc); Src.evals = SrcEigenvalues(1:numEigsSrc);
    Tar.evecs = TarLaplaceBasis(:,1:numEigsTar); Tar.evals = TarEigenvalues(1:numEigsTar);
    
    %% Descriptors
    fct_src = [];
    % fprintf('Computing the descriptors...\n');tic;
    fct_src = [fct_src, waveKernelSignature(SrcLaplaceBasis, SrcEigenvalues, Src.A, 200)];
    fct_src = [fct_src, waveKernelMap(SrcLaplaceBasis, SrcEigenvalues, Src.A, 200, landmarks(:,1))];
    
    fct_tar = [];
    fct_tar = [fct_tar, waveKernelSignature(TarLaplaceBasis, TarEigenvalues, Tar.A, 200)];
    fct_tar = [fct_tar, waveKernelMap(TarLaplaceBasis, TarEigenvalues, Tar.A, 200, landmarks(:,2))];
    
    % Subsample descriptors (for faster computation). More descriptors is
    % usually better, but can be slower.
    fct_src = fct_src(:,1:40:end);
    fct_tar = fct_tar(:,1:40:end);
    
    % fprintf('done computing descriptors (%d on source and %d on target)\n',size(fct_src,2),size(fct_tar,2)); toc;
    
    assert(size(fct_src,2)==size(fct_tar,2));
    
    % Normalization
    no = sqrt(diag(fct_src'*Src.A*fct_src))';
    fct_src = fct_src ./ repmat(no, [Src.nv,1]);
    fct_tar = fct_tar ./ repmat(no, [Tar.nv,1]);
    
    %    fprintf('Pre-computing the multiplication operators...');tic;
    %% Multiplication Operators
    numFct = size(fct_src,2);
    OpSrc = cell(numFct,1);
    OpTar = cell(numFct,1);
    for i = 1:numFct
        OpSrc{i} = Src.evecs'*Src.A*(repmat(fct_src(:,i), [1,numEigsSrc]).*Src.evecs);
        OpTar{i} = Tar.evecs'*Tar.A*(repmat(fct_tar(:,i), [1,numEigsTar]).*Tar.evecs);
    end
    
    Fct_src = Src.evecs'*Src.A*fct_src;
    Fct_tar = Tar.evecs'*Tar.A*fct_tar;
    %  fprintf('done\n');toc;
    
    %% Fmap Computation
    %fprintf('Optimizing the functional map...\n');tic;
    Dlb = (repmat(Src.evals, [1,numEigsTar]) - repmat(Tar.evals', [numEigsSrc,1])).^2;
    Dlb = Dlb/norm(Dlb, 'fro')^2;
    constFct = sign(Src.evecs(1,1)*Tar.evecs(1,1))*[sqrt(sum(Tar.area)/sum(Src.area)); zeros(numEigsTar-1,1)];
    
    Dlb2 = (repmat(Tar.evals, [1,numEigsSrc]) - repmat(Src.evals', [numEigsTar,1])).^2;
    Dlb2 = Dlb2/norm(Dlb2, 'fro')^2;
    constFct2 = sign(Tar.evecs(1,1)*Src.evecs(1,1))*[sqrt(sum(Src.area)/sum(Tar.area)); zeros(numEigsSrc-1,1)];
    
    a = 1e-1; % Descriptors preservation
    b = 0;    % Commutativity with descriptors
    c = 1e-3; % Commutativity with Laplacian
    funObj = @(F) deal( a*sum(sum((reshape(F, [numEigsTar,numEigsSrc])*Fct_src - Fct_tar).^2))/2 + b*sum(cell2mat(cellfun(@(X,Y) sum(sum((X*reshape(F, [numEigsTar,numEigsSrc]) - reshape(F, [numEigsTar,numEigsSrc])*Y).^2)), OpTar', OpSrc', 'UniformOutput', false)), 2)/2 + c*sum( (F.^2 .* Dlb(:))/2 ),...
        a*vec((reshape(F, [numEigsTar,numEigsSrc])*Fct_src - Fct_tar)*Fct_src') + b*sum(cell2mat(cellfun(@(X,Y) vec(X'*(X*reshape(F, [numEigsTar,numEigsSrc]) - reshape(F, [numEigsTar,numEigsSrc])*Y) - (X*reshape(F, [numEigsTar,numEigsSrc]) - reshape(F, [numEigsTar,numEigsSrc])*Y)*Y'), OpTar', OpSrc', 'UniformOutput', false)), 2) + c*F.*Dlb(:));
    funProj = @(F) [constFct; F(numEigsTar+1:end)];
    
    funObj2 = @(F) deal( a*sum(sum((reshape(F, [numEigsTar,numEigsSrc])*Fct_tar - Fct_src).^2))/2 + b*sum(cell2mat(cellfun(@(X,Y) sum(sum((X*reshape(F, [numEigsTar,numEigsSrc]) - reshape(F, [numEigsTar,numEigsSrc])*Y).^2)),  OpSrc', OpTar', 'UniformOutput', false)), 2)/2 + c*sum( (F.^2 .* Dlb2(:))/2 ),...
        a*vec((reshape(F, [numEigsTar,numEigsSrc])*Fct_tar - Fct_src)*Fct_tar') + b*sum(cell2mat(cellfun(@(X,Y) vec(X'*(X*reshape(F, [numEigsTar,numEigsSrc]) - reshape(F, [numEigsTar,numEigsSrc])*Y) - (X*reshape(F, [numEigsTar,numEigsSrc]) - reshape(F, [numEigsTar,numEigsSrc])*Y)*Y'), OpSrc', OpTar', 'UniformOutput', false)), 2) + c*F.*Dlb2(:));
    
    funProj2 = @(F) [constFct2; F(numEigsTar+1:end)];
    
    funProj3 = @(F) [funProj(F(1:end/2)); funProj2(F(end/2+1:end))];
    
    %%
    F_lb = zeros(numEigsTar*numEigsSrc, 1); F_lb(1) = constFct(1);
    F_lb2 = zeros(numEigsTar*numEigsSrc, 1); F_lb2(1) = constFct2(1);
    
    %%
    % Compute the optional functional map using a quasi-Newton method.
    options.verbose = 0;
    Finit = [F_lb; F_lb2];
    
    lb1 = diag(Src.evals(1:numEigsSrc));
    lb2 = diag(Tar.evals(1:numEigsTar));
    
    F = minConf_PQN(@(F) funObj4(F, epss1,  epss, numEigsSrc, numEigsTar, funObj, funObj2, lb1, lb2), Finit, funProj3, options);
    C1 = reshape(F(1:end/2),numEigsTar, numEigsSrc);
    C2 = reshape(F(end/2+1:end),numEigsTar, numEigsSrc);
    
    %%
    F_lb = C1;
    [F_lb2, ~] = icp_refine(Src.evecs, Tar.evecs, C1, 5);
    
    %% Evaluation
    % Compute the p2p map
    
    % fmap before ICP (for comparison)
    pF_lb = knnsearch((F_lb*Src.evecs')', Tar.evecs);
    % fmap after ICP
    pF_lb2 = knnsearch((F_lb2*Src.evecs')', Tar.evecs);
    
    map_Ad = pF_lb;
    map_AdICP =  pF_lb2;
    
    fps_src = dijkstra_fps(Tar, 300);
    % compute the errors
    fprintf('errors with adjoint:\n');

    errsA = dijkstra_pairs(Src,[pF_lb(fps_src) fps_src])/Src.sqrt_area;
    fprintf('Mean map error (without ICP): %f\n',mean(errsA));

    errsICPA = dijkstra_pairs(Src, [pF_lb2(fps_src) fps_src])/Src.sqrt_area;
    fprintf('Mean map error (with ICP): %f\n', mean(errsICPA));
    
    F = minConf_PQN(@(F) funObj4(F, 0,  0, numEigsSrc, numEigsTar, funObj, funObj2, lb1, lb2), Finit, funProj3, options);
    C1 = reshape(F(1:end/2),numEigsTar, numEigsSrc);
    
    %%
    F_lb = C1;
    [F_lb2, ~] = icp_refine(Src.evecs, Tar.evecs, F_lb, 5);
    
    %% Evaluation
    % Compute the p2p map
    
    % fmap before ICP (for comparison)
    pF_lb = knnsearch((F_lb*Src.evecs')', Tar.evecs);
    % fmap after ICP
    pF_lb2 = knnsearch((F_lb2*Src.evecs')', Tar.evecs);
    
    map_NA = pF_lb;
    map_NAICP = pF_lb2;
    
    % compute the errors
    fprintf('errors without adjoint:\n');
    
    errsNA = dijkstra_pairs(Src,[pF_lb(fps_src) fps_src])/Src.sqrt_area;
    fprintf('Mean map error (without ICP): %f\n',mean(errsNA));
    errsICPNA = dijkstra_pairs(Src, [pF_lb2(fps_src) fps_src])/Src.sqrt_area;
    fprintf('Mean map error (with ICP): %f\n', mean(errsICPNA));
    
    meanerrs = [meanerrs; mean(errsA) mean(errsICPA) mean(errsNA) mean(errsICPNA)];
    
    allerrsA = [allerrsA errsA];
    allerrsICPA = [allerrsICPA errsICPA];
    
    allerrsNA = [allerrsNA errsNA];
    allerrsICPNA = [allerrsICPNA errsICPNA];
    
    hold off;
    nr = length(reshape(allerrsA,[],1));
    plot(sort(reshape(allerrsA,[],1)),linspace(0,1,nr),'-g','LineWidth',2);
    hold on;
    plot(sort(reshape(allerrsICPA,[],1)),linspace(0,1,nr),'--g','LineWidth',2);
    
    plot(sort(reshape(allerrsNA,[],1)),linspace(0,1,nr),'-b','LineWidth',2);
    plot(sort(reshape(allerrsICPNA,[],1)),linspace(0,1,nr),'--b','LineWidth',2);
    axis([0 0.25 0 1]);
    pause(0.01);
end

FigHandle = figure('Position', [100, 100, 800, 600]);
hold on;
set(gca,'FontSize',20);
orangec = [1 0.7 0];
nr = length(reshape(allerrsCF,[],1));

title('TOSCA (76 pairs)','FontSize',24,'FontWeight','b');
plot(sort(reshape(allerrsICPA,[],1)),linspace(0,1,nr),'--','LineWidth',3,'Color',orangec);
plot(sort(reshape(allerrsA,[],1)),linspace(0,1,nr),'-','LineWidth',3,'Color',orangec);
plot(sort(reshape(allerrsICPCF,[],1)),linspace(0,1,nr),'--','LineWidth',3,'Color',[0 0.7 0.2]);
plot(sort(reshape(allerrsCF,[],1)),linspace(0,1,nr),'-','LineWidth',3,'Color',[0 0.7 0.2]);
plot(sort(reshape(allerrsICPNA,[],1)),linspace(0,1,nr),'--','LineWidth',3,'Color','b','MarkerSize',5);
plot(sort(reshape(allerrsNA,[],1)),linspace(0,1,nr),'-','LineWidth',3,'Color','b');
axis([0 0.25 0 1]);
xlabel('Geodesic Error','FontSize',24,'FontWeight','b');
ylabel('Fraction of Correspondences','FontSize',24,'FontWeight','b');
h_legend = legend('Adjoint Regularization + ICP','Adjoint Regularization',...
    '[ERGB] + ICP','[ERGB]','Regular Fmaps + ICP','Regular Fmaps',...
    'Location','southeast');
set(h_legend,'FontSize',22);
box on;

