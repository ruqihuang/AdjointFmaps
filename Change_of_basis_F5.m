clear all; close all;
addpath(genpath('data')); 
addpath(genpath('external')); 
addpath(genpath('utils')); 

S1 = read_off_shape('mesh000_simplified.off'); 
S2 = read_off_shape('mesh003_simplified.off'); 

S1 = compute_laplacian_basis(S1, 3); 
S2 = compute_laplacian_basis(S2, 752); 

samp = 17:15:752;
error1 = zeros(length(samp), 1); 
error2 = error1; 

for k = 1:length(samp)
    Phi2 = S2.evecs(:, 1:samp(k)); 
    C = Phi2\S1.evecs; 
    
    H = Phi2'*S1.A*Phi2; 
    [e,v] = eig(diag(S2.evals(1:samp(k))),H);
    e = e*diag(1./sqrt(diag(e'*H*e)));
    [v, order] = sort(diag(v));
    e = e(:,order);
    
    C2 = (Phi2*e)\S1.evecs; 
    
    [~, match1] = icp_refine(S1.evecs, Phi2, C, 0); 
    [~, match2] = icp_refine(S1.evecs, Phi2*e, C2, 0); 
    
    % compute overall error with respect to the two recovery methods
    error1(k) = sqrt(sum(sum((S1.Pts - S1.Pts(match1, :)).^2, 2))); 
    error2(k) = sqrt(sum(sum((S1.Pts - S1.Pts(match2, :)).^2, 2))); 
    
end

figure; 
plot(samp, error1, '-b','LineWidth',2); 
hold on; 
plot(samp, error2, '-r', 'LineWidth', 2); 
hold on; 
scatter(samp, error1, 50, 'filled'); 
hold on; 
scatter(samp, error2, 50, 'filled'); 
xlabel('Number of eigenfunctions used in the target shape'); 
ylabel('Recovery error'); 
legend('Without a change of basis', 'With a change of basis', 'Location', 'east'); 
    
