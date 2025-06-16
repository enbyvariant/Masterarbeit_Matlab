function [A, b, H, c, c_0, xl, xu, C, cl, cu] = get_input()
    
    prob = cutest_setup();
    n = prob.n;
    m = prob.m;
    r = 0;
    
    % compute correct dimensions for m and r
    for i = 1:m
        if prob.cl(i) || prob.cu(i)
            r = r + 1;
            m = m - 1;
        end
    end
    
    % variables directly available from cutest
    const = -cutest_cons(zeros(n,1));
    c_0 = cutest_obj(zeros(n,1));
    xl = sparse(prob.bl);
    xu = sparse(prob.bu);
    H = cutest_hess(prob.x, prob.v);

    M = zeros(m,n);
    c = zeros(n,1);
    for i = 1:n
        M(1:m,i) = cutest_cons([zeros(i-1,1); 1; zeros(n-i,1)])+const;
        c(i) = cutest_obj([zeros(i-1,1); 1; zeros(n-i,1)])- 1/2*H(i,i)-c_0;
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

    for j = 1:m
        % constraint of the form Cx <= cu
        if prob.cl(j) < - 10^3
            i = i + 1;
            C(i,1:n) = M(j,1:n);
            cu(i) = const(j);
        end

        % constraint of the form cl <= Cx
        if prob.cu > 10^3
            i = i + 1;
            C(i,1:n) = M(j,1:n);
            cl(i) = const(j);
        end

        % constraint of the form Ax = b
        if ~prob.cl(j) && ~prob.cu(j)
            k = k + 1;
            A(k,1:n) = M(j,1:n);
            b(k) = const(j);
        end
    end
    
    H = sparse(H);
    c = sparse(c);

    C = sparse(C);
    cl = sparse(cl);
    cu = sparse(cu);

    A = sparse(A);
    b = sparse(b);

    

    if any(prob.cl)
        fprintf("error in input: c_l =");
        disp(prob.cl);
    end
    if any(prob.cu)
        disp("error in input: c_u =")
        disp(prob.cu);
    end

    cutest_terminate;
end