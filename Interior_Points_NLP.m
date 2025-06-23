function [outputArg1,outputArg2] = Interior_Points_NLP(x,inputArg2)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    p = cutest_setup;
    n = p.n;
    m = prob.m;
    r = 0;
    
    % compute correct dimensions for m and r
    for i = 1:m
        if prob.cl(i) || prob.cu(i)
        r = r + 1;
        m = m - 1;
        end
    end

    for j = 1:m+r
        % constraint of the form c(x) <= cu
        if prob.cl(j) < - 10^3
            i = i + 1;
            C(i,1:n) = M(j,1:n);
            cu(i) = const(j);
            cl(i) = - 10^5;
        end

        % constraint of the form cl <= c(x)
        if prob.cu(j) > 10^3
            i = i + 1;
            C(i,1:n) = M(j,1:n);
            cl(i) = const(j);
            cu(i) = 10^5;
        end

        % constraint of the form c(x) = 0
        if ~prob.cl(j) && ~prob.cu(j)
            k = k + 1;
            A(k,1:n) = M(j,1:n);
            b(k) = const(j);
        end
    end


    while 1
        grad = cutest_grad(x);
    
    end

outputArg1 = inputArg1;
outputArg2 = inputArg2;
end