function [nlp, dim] = get_input()
    
    prob = cutest_setup();
    n = prob.n;
    m = prob.m;
    r = 0;
    
    % compute correct dimensions for m and r
    for i = 1:prob.m
        if prob.cl(i) || prob.cu(i)
            r = r + 1;
            m = m - 1;
        end
    end
    dim = struct('n', n, 'm', m, 'r', r);
    
    % variables directly available from cutest
    const = -cutest_cons(zeros(n,1));
    nlp.c_0 = cutest_obj(zeros(n,1));
    nlp.xl = sparse(prob.bl);
    nlp.xu = sparse(prob.bu);
    H = cutest_hess(prob.x, prob.v);
    
    % compute all constraint multipliers and the linear part of obj func
    M = zeros(m+r,n);
    c = zeros(n,1);
    for i = 1:n
        M(1:m+r,i) = cutest_cons([zeros(i-1,1); 1; zeros(n-i,1)])+const;
        c(i) = cutest_obj([zeros(i-1,1); 1; zeros(n-i,1)])- 1/2*H(i,i)-nlp.c_0;
    end

    % variables for nontrivial inequalities
    i = 0;
    C = zeros(r,n);
    cl = zeros(r,1);
    cu = zeros(r,1);
    % variables for equalities
    k = 0;
    A = zeros(m,n);
    b = zeros(m,1);

    for j = 1:m+r
        % constraint of the form Cx <= cu
        if prob.cl(j) < - 10^3
            i = i + 1;
            C(i,1:n) = M(j,1:n);
            cu(i) = const(j);
            cl(i) = - 10^5;
        end

        % constraint of the form cl <= Cx
        if prob.cu(j) > 10^3
            i = i + 1;
            C(i,1:n) = M(j,1:n);
            cl(i) = const(j);
            cu(i) = 10^5;
        end

        % constraint of the form Ax = b
        if ~prob.cl(j) && ~prob.cu(j)
            k = k + 1;
            A(k,1:n) = M(j,1:n);
            b(k) = const(j);
        end
    end
    
    nlp.H = sparse(H);
    nlp.c = sparse(c);

    nlp.C = sparse(C);
    nlp.cl = sparse(cl);
    nlp.cu = sparse(cu);

    nlp.A = sparse(A);
    nlp.b = sparse(b);

    cutest_terminate;
end