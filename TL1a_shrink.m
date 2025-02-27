
function [Y,rankY] = TL1Norm(M,mu,m,n)
    a=1;
    [U, S, V] = svd(M);
    S_shrunk = diag(shrinkTL1(diag(S),mu,a));
    S_new = zeros(m,n); S_new(1:min(m,n),1:min(m,n)) = S_shrunk;
    Y = U * S_new * V';
    rankY = nnz(S_new);
    
%=============================================================================

function v = shrinkTL1(s,lambda,a)
    phi = acos(1-(0.5*27*lambda*a*(a+1))./(a+abs(s)).^3);
    v = sign(s).*(2/3 * (a+abs(s)).* cos(phi/3) -2*a/3+abs(s)/3).*(abs(s)>lambda);
