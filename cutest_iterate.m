function [grad, H, cons_e, cons_i, cons_e_grad, cons_i_grad] = cutest_iterate(x,y,sl,su,tl,tu,wl,wu,zl,zu,n,m)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
p = cutest_setup;
j = 0;
k = 0;
i = 0;
equ = zeros(m);
for i = 1:p.m
    if p.cl(i) == 0 && p.cu(i) == 0
        j = j+1;
        equ(j) = i;
    else
        k = k+1;
        inequ(k) = i;
    end


end


grad = cutest_grad(x);
H = cutest_hess(x,lambda);
outputArg2 = inputArg2;
end