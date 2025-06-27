function [grad, H, cons_e, cons_i, cons_e_grad, cons_i_grad, obj] = cutest_iterate(x,y,zl,zu,equ,inequ)
% Gives out the values for the current iterate
%  Using the cutest environment
    n = size(x,1);
    m = size(equ,1);
    r = size(inequ,1);

    p = cutest_setup;
    
    [cons, C] = cutest_cons(x);
    cons_e = cons(equ);
    cons_i = cons(inequ);
    cons_e_grad = C(equ,equ);
    cons_i_grad = C(inequ,inequ);

    lambda = zeros(m+r,1);
    lambda(equ) = y;
    lambda(inequ) = zl - zu;

    [L, grad] = cutest_lag(x,lambda);
    H = cutest_hess(x,lambda);
    
    obj = cutest_obj(x);

    cutest_terminate;
end