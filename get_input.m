function [nlp, dim, f] = get_input()
    
    prob = cutest_setup();
    dim.n = prob.n;
    dim.m = prob.m;
    dim.r = 0;
    
    f = 0;
    % if (dim.n < 1000) || (dim.n > 7000 || dim.m > 7000)
    %     f = 1;
    %     nlp = [];
    %     fprintf(1, 'n = %i\n', dim.n);
    %     fprintf(1, 'm = %i\n', dim.m);
    %     cutest_terminate;
    %     return
    % end
    
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
    H = cutest_sphess(prob.x, prob.v);

    c = zeros(dim.n,1);
    nz = 0;
    nz_i = zeros(dim.n,1);
    for k = 1:dim.n
        name = genvarname(string(k));
        cons_column.(name) = cutest_cons([zeros(k-1,1); 1; zeros(dim.n-k,1)]);
        cons_column.(name) = cons_column.(name) + const;
        nz_i(k) = nnz(cons_column.(name));
        nz = nz + nz_i(k);
    end
    i = zeros(nz,1);
    j = zeros(nz,1);
    v = zeros(nz,1);
    l = 1;
    for k = 1:dim.n
        name = genvarname(string(k));
        [ic,jc,vc] = find(cons_column.(name));
        i(l:l + nz_i(k) -1,1) = ic;
        j(l:l + nz_i(k) -1,1) = ones(size(jc)) * k;
        v(l:l + nz_i(k) -1,1) = vc;
        l = l + nz_i(k);
    end


    M = sparse(i,j,v, dim.m + dim.r,dim.n);

    for i = 1:dim.n
        c(i) = cutest_obj([zeros(i-1,1); 1; zeros(dim.n-i,1)])- 1/2*H(i,i)-c_0;
    end

    % variables for nontrivial inequalities
    n_1 = 0;
    cl = zeros(dim.r,1);
    cu = zeros(dim.r,1);
    % variables for equalities
    n_2 = 0;
    b = zeros(dim.m,1);
    [i,j,v] = find(M);
    nz = size(v);
    Ci = zeros(nz);
    Cj = zeros(nz);
    Cv = zeros(nz);
    Ai = zeros(nz);
    Aj = zeros(nz);
    Av = zeros(nz);
    nz_A = 0;
    nz_C = 0;

    for k = 1:nz(1)
        if prob.cl(i(k)) || prob.cu(i(k))
            Ci(k) = i(k) - dim.m;
            Cj(k) = j(k);
            Cv(k) = v(k);
            nz_C = nz_C + 1;
        else
            Ai(k) = i(k);
            Aj(k) = j(k);
            Av(k) = v(k);
            nz_A = nz_A + 1;
        end
    end
    l_1 = 1;
    l_2 = 1;
    if dim.r == 0
        Ci = [];
        Cj = [];
        Cv = [];
    else
        while l_1 < length(Ci) + 1
            if Ci(l_1) == 0
                Ci(l_1) = [];
                Cj(l_1) = [];
                Cv(l_1) = [];
            else
                l_1 = l_1 + 1;
            end
        end
    end
    if dim.m == 0
        Ai = [];
        Aj = [];
        Av = [];
    end
    while l_2 < length(Ai) + 1
        if Ai(l_2) == 0
            Ai(l_2) = [];
            Aj(l_2) = [];
            Av(l_2) = [];
        else
            l_2 = l_2 + 1;
        end
    end
    C = sparse(Ci,Cj,Cv,dim.r,dim.n);
    A = sparse(Ai,Aj,Av,dim.m,dim.n);
             

    for j = 1:dim.m+dim.r
        % constraint of the form Cx <= cu
        if prob.cl(j) < - 10^7
            n_1 = n_1 + 1;
            cu(n_1) = const(j);
            cl(n_1) = -10^7;
        end

        % constraint of the form cl <= Cx
        if prob.cu(j) > 10^7
            n_1 = n_1 + 1;
            cl(n_1) = const(j);
            cu(n_1) = 10^7;
        end

        % constraint of the form Ax = b
        if ~prob.cl(j) && ~prob.cu(j)
            n_2 = n_2 + 1;
            b(n_2) = const(j);
        end
    end

    nlp.index_xl = find(xl > -10^7);
    nlp.index_xu = find(xu < 10^7);
    nlp.index_cl = find(cl > -10^7);
    nlp.index_cu = find(cu < 10^7);
    
    nlp.H = H;
    nlp.c = sparse(c);
    nlp.c_0 = c_0;

    nlp.C = C;
    nlp.cl = sparse(cl);
    nlp.cu = sparse(cu);

    nlp.xl = xl;
    nlp.xu = xu;

    nlp.A = A;
    nlp.b = sparse(b);

    cutest_terminate;
end