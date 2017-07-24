function S = compute_laplacian_basis(S1, numEigs)

	S.surface = S1.surface;

    X = [S1.surface.X S1.surface.Y S1.surface.Z];
    T = S1.surface.TRIV; 

%     % compute face normal
%     
%     S.Nf = cross(X(T(:,1),:) - X(T(:,2),:), X(T(:,1),:) - X(T(:,3),:));
% 	S.Nf = S.Nf./repmat(sqrt(sum(S.Nf.^2, 2)), [1, 3]);

    S.A = diag(vertexAreas(X, S1.surface.TRIV));
    S.W = cotWeights(X, S1.surface.TRIV);
	
    nv = size(S.A,1);

	% compute laplacian eigenbasis.
	try
	    [S.evecs, S.evals] = eigs(S.W, S.A, numEigs, 1e-6);
	catch
	    % In case of trouble make the laplacian definite
	    [S.evecs, S.evals] = eigs(S.W - 1e-8*speye(nv), S.A, numEigs, 'sm');
	end

    [S.evals, order] = sort(diag(S.evals), 'ascend');
    S.evecs = S.evecs(:,order);

	S.nv = nv; 
    S.Pts = X; 
	S.area = diag(S.A); 
    S.sqrt_area = sqrt(sum(S.area));  
    S = normals(S); 
    
    % to be compatible
    S.X = S.Pts; 
    S.T = S.surface.TRIV;
end
