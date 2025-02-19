
function v = shrinkTL1(s,lambda,a)
    phi = acos(1-(0.5*27*lambda*a*(a+1))./(a+abs(s)).^3);
    v = sign(s).*(2/3 * (a+abs(s)).* cos(phi/3) -2*a/3+abs(s)/3).*(abs(s)>lambda);