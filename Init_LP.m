function [x, sl, su, tl, tu, y, wl, wu, zl, zu, bound_xl, bound_xu, bound_cl, bound_cu] = Init_LP(nlp,dim)
%
% 
    % Initial values
    n = dim.n;
    m = dim.m;
    r = dim.r;
    en = ones(n,1);
    er = ones(r,1);
    
    index_l = ones(1,n);
    index_u = ones(1,n);
    wl = ones(n,1);
    wu = ones(n,1);
    for i = 1:n
        if nlp.xl(i) < -10^3
            nlp.xl(i) = 0;
            wl(i) = 0; 
            index_l(i) = 0;
        end
        if nlp.xu(i) > 10^3
            nlp.xu(i) = 0;
            wu(i) = 0;
            index_u(i) = 0;
        end
    end
    index_l = find(index_l);
    index_u = find(index_u);
    bound_xl = sparse(index_l, index_l, ones(size(index_l)), n, n);
    bound_xu = sparse(index_u, index_u, ones(size(index_u)), n, n);

    
    index_cl = ones(1,r);
    index_cu = ones(1,r);
    zl = ones(r,1);
    zu = ones(r,1);
    for i = 1:r
        if nlp.cl(i) < -10^3
            nlp.cl(i) = 0;
            zl(i) = 0;
            index_cl(i) = 0;
        end
        if nlp.cu(i) > 10^3
            nlp.cu(i) = 0;
            zu(i) = 0;
            index_cu(i) = 0;
        end
    end
    index_cl = find(index_cl);
    index_cu = find(index_cu);
    bound_cl = sparse(index_cl, index_cl, ones(size(index_cl)), r, r);
    bound_cu = sparse(index_cu, index_cu, ones(size(index_cu)), r, r);

    I_n = eye(n);
    I_r = eye(r);
    n_l = size(nlp.index_xl,1);
    n_u = size(nlp.index_xu,1);
    r_l = size(nlp.index_cl,1);
    r_u = size(nlp.index_cu,1);

    M = [nlp.A zeros(m, n_l + n_u + r_l + r_u);
        I_n(nlp.index_xl,:)  -I_n(nlp.index_xl,nlp.index_xl)       zeros(n_l,n_u)                zeros(n_l,r_l) zeros(n_l,r_u);
        -I_n(nlp.index_xu,:)     zeros(n_u,n_l)                -I_n(nlp.index_xu,nlp.index_xu)   zeros(n_u,r_l) zeros(n_u,r_u);
        nlp.C(nlp.index_cl,:)    zeros(r_l,n_l)                       zeros(r_l,n_u)             -I_r(nlp.index_cl,nlp.index_cl) zeros(r_l,r_u);
        -nlp.C(nlp.index_cu,:)   zeros(r_u,n_l)                       zeros(r_u,n_u)            zeros(r_u,r_l) -I_r(nlp.index_cu,nlp.index_cu)
        ];

    % Compute starting point
    vector = [nlp.b; nlp.xl(nlp.index_xl); nlp.xu(nlp.index_xu); nlp.cl(nlp.index_cl); nlp.cu(nlp.index_cu)];

    a_prim = M' * (M*M'\vector);
    x = a_prim(1:n);
    sl = bound_xl*(x - nlp.xl);
    su = bound_xu*(-x + nlp.xu);
    tl = bound_cl*(nlp.C*x - nlp.cl);
    tu = bound_cu*(-nlp.C*x + nlp.cu);
    
    N = [nlp.A' eye(n) nlp.C'];
    a_dual = N' * (N*N'\nlp.c);
    y = a_dual(1:m);
    wl = a_dual(m+1:m+n);
    wu = zeros(n,1);

    for i =1:r
        if ismember(i,nlp.index_cl)
            zl(i) = a_dual(m+n+i);
            zu(i) = 0;
        else
            zu(i) = a_dual(m+n+i);
            zl(i) = 0;
        end
    end
    
    delta_pri = 3/2*max([-sl;-su;-tl;-tu;0]);
    delta_dual = 3/2*max([-wl;-wu;-zl;-zu; 0]);
    
    sl = bound_xl*(sl + delta_pri * en);
    su = bound_xu*(su + delta_pri * en);
    tl = bound_cl*(tl + delta_pri * er);
    tu = bound_cu*(tu + delta_pri * er);
    wl = bound_xl*(wl + delta_dual * en);
    wu = bound_xu*(wu + delta_dual * en);
    zl = bound_cl*(zl + delta_dual * er);
    zu = bound_cu*(zu + delta_dual * er);

    delta_pri = 1/2*(sl' * wl + su'*wu + tl'*zl + tu'*zu)/(en' * wl + en'*wu + er'*zl + er'*zu);
    delta_dual = 1/2*(sl' * wl + su'*wu + tl'*zl + tu'*zu)/(en' * sl + en'*su + er'*tl + er'*tu);

    sl = bound_xl*(sl + delta_pri * en);
    su = bound_xu*(su + delta_pri * en);
    tl = bound_cl*(tl + delta_pri * er);
    tu = bound_cu*(tu + delta_pri * er);
    wl = bound_xl*(wl + delta_dual * en);
    wu = bound_xu*(wu + delta_dual * en);
    zl = bound_cl*(zl + delta_dual * er);
    zu = bound_cu*(zu + delta_dual * er);
    

    % Failsave, if some linear equation system is impossible to solve
    if any(isnan(x))
        x = ones(n,1);
    end
    if any(isnan(sl))
        sl = bound_xl*ones(n,1);
    end
    if any(isnan(su))    
        su = bound_xu*ones(n,1);
    end
    if any(isnan(tl))
        tl = bound_cl*ones(r,1);
    end
    if any(isnan(tu))
        tu = bound_cu*ones(r,1);
    end
    if any(isnan(wl))
        wl = bound_xl*ones(n,1);
    end

    if any(isnan(wu))
        wu = bound_xu*ones(n,1);
    end
    if any(isnan(zl))
        zl = bound_cl*ones(r,1);
    end
    if any(isnan(zu))
        zu = bound_cu*ones(r,1);
    end
    if any(isnan(y))
        y = ones(m,1);
    end

end