function [nlp, obj] = cutest_iterate(it,nlp, dim)
% Gives out the values for the current iterate
%  Using the cutest environment

    p = cutest_setup;
    
    [cons, Jacobian] = cutest_cons(it.x);
    nlp.cons_e = cons(nlp.equ);
    nlp.cons_i = cons(nlp.inequ);
    nlp.A = Jacobian(nlp.equ,nlp.equ);
    nlp.C = Jacobian(nlp.inequ,nlp.inequ);

    lambda = zeros(dim.m+dim.r,1);
    lambda(nlp.equ) = it.y;
    lambda(nlp.inequ) = it.zl - it.zu;

    [nlp.L, nlp.grad] = cutest_lag(it.x,lambda);
    nlp.H = cutest_hess(it.x,lambda);
    
    obj = cutest_obj(it.x);

    cutest_terminate;
end