function [a,b] = funObj4(F,eps1, eps2, numEigsTar, numEigsSrc,funObj,funObj2,lb1,lb2)
[a1,b1] = funObj(F(1:end/2));
[a2,b2] = funObj2(F(end/2+1:end));

C1 = reshape(F(1:end/2),numEigsTar, numEigsSrc);
C2 = reshape(F(end/2+1:end),numEigsTar, numEigsSrc);

a = a1 + a2 + eps1*norm(C1-C2','fro').^2/2 + eps2*norm(lb2*C1-C2'*lb1,'fro').^2/2;

b = [b1 + eps1*reshape(C1-C2',[],1) + eps2*reshape(lb2*(lb2*C1-C2'*lb1),[],1);...
    b2 + eps1*reshape(C2-C1',[],1) + eps2*reshape(lb1*(lb1*C2-C1'*lb2),[],1)];


% b = [zeros(size(b2));b2];