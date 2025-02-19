
function [Y,rankY] = TL1Norm(M,mu,m,n)
    a=100;
    [U, S, V] = svd(M);
    S_shrunk = diag(shrinkTL1(diag(S),mu,a));
    S_new = zeros(m,n); S_new(1:min(m,n),1:min(m,n)) = S_shrunk;
    Y = U * S_new * V';
    rankY = nnz(S_new);