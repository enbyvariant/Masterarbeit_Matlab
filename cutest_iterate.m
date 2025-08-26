function [nlp] = cutest_iterate(it, nlp, dim,p)
% Gives out the values for the current iterate
%  Using the cutest environment
    m = dim.m;
    r = dim.r;
    n = dim.n;
    
    %disp(it.x);
    [cons, J] = cutest_cons(it.x);
    nlp.c_e = cons(nlp.equ);
    nlp.c_i = cons(nlp.inequ);
    nlp.A = J(nlp.equ,:);
    if isempty(nlp.A)
        nlp.A = zeros(m,n);
    end
    nlp.C = J(nlp.inequ,:);

    lambda = zeros(m+r,1);
    lambda(nlp.equ) = it.y;
    lambda(nlp.inequ) = it.zl - it.zu;

    [nlp.L, nlp.l_grad] = cutest_lag(it.x,lambda);
    nlp.H = cutest_hess(it.x,lambda);
    
    [nlp.obj, nlp.grad] = cutest_obj(it.x);

end