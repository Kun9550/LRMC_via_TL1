
function [X,Y,W] = TL1(A,options)
    sigma = 5e-3; 
    lamfactor = options.lamfactor;
    if isfield(options,'sigma'); sigma = options.sigma; end    
    
    Fnorm = @(X) mexFnorm(X); 
    tstart = clock;
    maxiter = 200; 
    stoptol = 1e-6;
    [m,n] = size(A);
    nzidx = find(abs(A) > 0); 
    Aobs = A(nzidx); 
    lamtarget = lamfactor*norm(Aobs);    
    Y = zeros(m,n);
    W = zeros(m,n);    
    breakyes = 0;
    tau = 1.618;     
    lam = 100*lamtarget; 
    
    for iter = 1:maxiter
        Yold = Y;        
        invsigma = 1/sigma;  
        lam = max(lamtarget,lam*0.95);
        Xinput = Y-invsigma*W; 
       [X,rankX] = TL1Norm(Xinput,lam*invsigma,m,n);
        Y = X+invsigma*W;
        Y(nzidx) = (1/(1+sigma))*(Aobs+sigma*Y(nzidx));
        W = W + tau*sigma*(X-Y);  
        
        normX = max(1,norm(X,'fro'));
        primfeas = norm(Y-X,'fro')/normX;
        dualfeas = Fnorm(Yold-Y) + invsigma*norm(X(nzidx)-Y(nzidx));
        if (max(primfeas,dualfeas) < stoptol) 
           if (lam==lamtarget)
              breakyes=1;
           else
              lam = lamtarget; 
           end
        end        
        if (rem(iter,20)==1) || (breakyes)
           ttime = etime(clock,tstart);                 
           fprintf(' %3.0f %3.2e %3.2e| %3.2e %3.0f| %5.2f\n',...
           iter,primfeas,dualfeas,sigma,rankX,ttime);  
        end
        if (rem(iter,10)==0) 
           if (primfeas < 0.5*dualfeas); 
              sigma = 1.2*sigma;
           elseif (primfeas > 2*dualfeas)
              sigma = 0.8*sigma; 
           end
        end
        if (breakyes); break; end
    end

%%************************************************************
%%************************************************************

function [Y,rankY] = TL1Norm(M,mu,m,n)
    a=100;
    [U, S, V] = svd(M);
    S_shrunk = diag(shrinkTL1(diag(S),mu,a));
    S_new = zeros(m,n); S_new(1:min(m,n),1:min(m,n)) = S_shrunk;
    Y = U * S_new * V';
    rankY = nnz(S_new);

%%************************************************************
%%************************************************************

function v = shrinkTL1(s,lambda,a)
    phi = acos(1-(0.5*27*lambda*a*(a+1))./(a+abs(s)).^3);
    v = sign(s).*(2/3 * (a+abs(s)).* cos(phi/3) -2*a/3+abs(s)/3).*(abs(s)>lambda);

%%************************************************************