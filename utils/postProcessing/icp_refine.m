function [C, matches] = icp_refine(L1, L2, C, nk)
    n1 = size(C,2);
    n2 = size(C,1);
        
    Vs = L1(:,1:n1);
    Vt = L2(:,1:n2);
    
    if(nk < 1)
        matches = knnsearch((C*Vs')',Vt);
    else
        for k=1:nk
            matches = knnsearch((C*Vs')',Vt);
            W = L2\L1(matches,:);
            if(k<nk-1)
                [s,~,d] = svd(W(1:n2,1:n1));
                C = s*eye(n2,n1)*d';
            end
        end
    end
end