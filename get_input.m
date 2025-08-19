function [nlp, dim] = get_input()
    
    prob = cutest_setup();
    dim.n = prob.n;
    dim.m = prob.m;
    dim.r = 0;
    
    % compute correct dimensions for m and r
    for i = 1:prob.m
        if prob.cl(i) || prob.cu(i)
            dim.r = dim.r + 1;
            dim.m = dim.m - 1;
        end
    end    

    % variables directly available from cutest
    const = -cutest_cons(zeros(dim.n,1));
    c_0 = cutest_obj(zeros(dim.n,1));
    xl = sparse(prob.bl);
    xu = sparse(prob.bu);
    H = cutest_hess(prob.x, prob.v);

    M = zeros(dim.m+dim.r,dim.n);
    c = zeros(dim.n,1);
    for i = 1:dim.n
        M(1:dim.m+dim.r,i) = cutest_cons([zeros(i-1,1); 1; zeros(dim.n-i,1)])+const;
        c(i) = cutest_obj([zeros(i-1,1); 1; zeros(dim.n-i,1)])- 1/2*H(i,i)-c_0;
    end

    % variables for nontrivial inequalities
    i = 0;
    C = zeros(dim.r,dim.n);
    cl = zeros(dim.r,1);
    cu = zeros(dim.r,1);
    % variables for equalities
    k = 0;
    A = zeros(dim.m,dim.n);
    b = zeros(dim.m,1);

    for j = 1:dim.m+dim.r
        % constraint of the form Cx <= cu
        if prob.cl(j) < - 10^3
            i = i + 1;
            C(i,1:dim.n) = M(j,1:dim.n);
            cu(i) = const(j);
            cl(i) = - 10^5;
        end

        % constraint of the form cl <= Cx
        if prob.cu(j) > 10^3
            i = i + 1;
            C(i,1:dim.n) = M(j,1:dim.n);
            cl(i) = const(j);
            cu(i) = 10^5;
        end

        % constraint of the form Ax = b
        if ~prob.cl(j) && ~prob.cu(j)
            k = k + 1;
            A(k,1:dim.n) = M(j,1:dim.n);
            b(k) = const(j);
        end
    end

    nlp.index_xl = find(xl > -10^3);
    nlp.index_xu = find(xu < 10^3);
    nlp.index_cl = find(cl > -10^3);
    nlp.index_cu = find(cu < 10^3);
    
    nlp.H = sparse(H);
    nlp.c = sparse(c);
    nlp.c_0 = c_0;

    nlp.C = sparse(C);
    nlp.cl = sparse(cl);
    nlp.cu = sparse(cu);

    nlp.xl = xl;
    nlp.xu = xu;

    nlp.A = sparse(A);
    nlp.b = sparse(b);

    cutest_terminate;
end