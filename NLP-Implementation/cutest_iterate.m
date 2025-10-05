function [nlp] = cutest_iterate(it, nlp, dim,p)
% Gives out the values for the current iterate
%  Using the cutest environment
    m = dim.m;
    r = dim.r;
    n = dim.n;
    
    [cons, J] = cutest_scons(it.x);
    nlp.c_e = cons(nlp.equ);
    nlp.c_i = cons(nlp.inequ);
    nlp.A = J(nlp.equ,:);
    if isempty(nlp.A)
        nlp.A = sparse(m,n);
    end
    nlp.C = J(nlp.inequ,:);

    lambda = zeros(m+r,1);
    lambda(nlp.equ) = it.y;
    lambda(nlp.inequ) = it.zl - it.zu;

    [nlp.L, nlp.l_grad] = cutest_lag(it.x,lambda);
    nlp.H = cutest_sphess(it.x,lambda);
    
    [nlp.obj, nlp.grad] = cutest_obj(it.x);
    [i,j,v] = find(nlp.A);
    for k = 1:length(v)
        if v(k) > -1*10^(-100) && v(k) < 1*10^(-100)
            v(k) = 0;
        end
    end
    nlp.A = sparse(i,j,v,m,n);

end
