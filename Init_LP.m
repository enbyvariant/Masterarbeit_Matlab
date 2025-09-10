function [it] = Init_LP(nlp,dim)
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
        if nlp.xl(i) < -10^7
            nlp.xl(i) = 0;
            wl(i) = 0; 
            index_l(i) = 0;
        end
        if nlp.xu(i) > 10^7
            nlp.xu(i) = 0;
            wu(i) = 0;
            index_u(i) = 0;
        end
    end
    index_l = find(index_l);
    index_u = find(index_u);
    it.bound_xl = sparse(index_l, index_l, ones(size(index_l)), n, n);
    it.bound_xu = sparse(index_u, index_u, ones(size(index_u)), n, n);

    
    index_cl = ones(1,r);
    index_cu = ones(1,r);
    zl = ones(r,1);
    zu = ones(r,1);
    for i = 1:r
        if nlp.cl(i) < -10^5
            nlp.cl(i) = 0;
            zl(i) = 0;
            index_cl(i) = 0;
        end
        if nlp.cu(i) > 10^5
            nlp.cu(i) = 0;
            zu(i) = 0;
            index_cu(i) = 0;
        end
    end
    index_cl = find(index_cl);
    index_cu = find(index_cu);
    it.bound_cl = sparse(index_cl, index_cl, ones(size(index_cl)), r, r);
    it.bound_cu = sparse(index_cu, index_cu, ones(size(index_cu)), r, r);
    
    n_i = 1:n;
    r_i = 1:r;
    I_n = sparse(n_i, n_i,ones(n,1));
    I_r = sparse(r_i,r_i,ones(r,1));
    n_l = size(nlp.index_xl,1);
    n_u = size(nlp.index_xu,1);
    r_l = size(nlp.index_cl,1);
    r_u = size(nlp.index_cu,1);

    M = [nlp.A sparse(m, n_l + n_u + r_l + r_u);
        I_n(nlp.index_xl,:)  -I_n(nlp.index_xl,nlp.index_xl)       sparse(n_l,n_u)                sparse(n_l,r_l) sparse(n_l,r_u);
        -I_n(nlp.index_xu,:)     sparse(n_u,n_l)                -I_n(nlp.index_xu,nlp.index_xu)   sparse(n_u,r_l) sparse(n_u,r_u);
        nlp.C(nlp.index_cl,:)    sparse(r_l,n_l)                       sparse(r_l,n_u)             -I_r(nlp.index_cl,nlp.index_cl) sparse(r_l,r_u);
        -nlp.C(nlp.index_cu,:)   sparse(r_u,n_l)                       sparse(r_u,n_u)            sparse(r_u,r_l) -I_r(nlp.index_cu,nlp.index_cu)
        ];

    % Compute starting point
    vector = [nlp.b; nlp.xl(nlp.index_xl); nlp.xu(nlp.index_xu); nlp.cl(nlp.index_cl); nlp.cu(nlp.index_cu)];
    sol = (M*M'\vector);
    a_prim = M' * sol;
    x = a_prim(1:n);
    sl = it.bound_xl*(x - nlp.xl);
    su = it.bound_xu*(-x + nlp.xu);
    tl = it.bound_cl*(nlp.C*x - nlp.cl);
    tu = it.bound_cu*(-nlp.C*x + nlp.cu);
    
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
    
    delta_pri = 3/2*max([-sl;-su;-tl;-tu;0.01]);
    delta_dual = 3/2*max([-wl;-wu;-zl;-zu; 0.01]);
    
    sl = it.bound_xl*(sl + delta_pri * en);
    su = it.bound_xu*(su + delta_pri * en);
    tl = it.bound_cl*(tl + delta_pri * er);
    tu = it.bound_cu*(tu + delta_pri * er);
    wl = it.bound_xl*(wl + delta_dual * en);
    wu = it.bound_xu*(wu + delta_dual * en);
    zl = it.bound_cl*(zl + delta_dual * er);
    zu = it.bound_cu*(zu + delta_dual * er);

    delta_pri = 1/2*(sl' * wl + su'*wu + tl'*zl + tu'*zu)/(en' * wl + en'*wu + er'*zl + er'*zu);
    delta_dual = 1/2*(sl' * wl + su'*wu + tl'*zl + tu'*zu)/(en' * sl + en'*su + er'*tl + er'*tu);

    sl = it.bound_xl*(sl + delta_pri * en);
    su = it.bound_xu*(su + delta_pri * en);
    tl = it.bound_cl*(tl + delta_pri * er);
    tu = it.bound_cu*(tu + delta_pri * er);
    wl = it.bound_xl*(wl + delta_dual * en);
    wu = it.bound_xu*(wu + delta_dual * en);
    zl = it.bound_cl*(zl + delta_dual * er);
    zu = it.bound_cu*(zu + delta_dual * er);
    

    % Failsave, if some linear equation system is impossible to solve
    if any(isnan(x))
        x = ones(n,1);
    end
    if any(isnan(sl))
        sl = it.bound_xl*ones(n,1);
    end
    if any(isnan(su))    
        su = it.bound_xu*ones(n,1);
    end
    if any(isnan(tl))
        tl = it.bound_cl*ones(r,1);
    end
    if any(isnan(tu))
        tu = it.bound_cu*ones(r,1);
    end
    if any(isnan(wl))
        wl = it.bound_xl*ones(n,1);
    end

    if any(isnan(wu))
        wu = it.bound_xu*ones(n,1);
    end
    if any(isnan(zl))
        zl = it.bound_cl*ones(r,1);
    end
    if any(isnan(zu))
        zu = it.bound_cu*ones(r,1);
    end
    if any(isnan(y))
        y = ones(m,1);
    end

    it.x = x;
    it.sl = sl;
    it.su = su;
    it.tl = tl;
    it.tu = tu;
    it.y = y;
    it.wl = wl;
    it.wu = wu;
    it.zl = zl;
    it.zu = zu;
    it.mu = (it.sl'*it.wl + it.su'*it.wu + it.tl'*it.zl + it.tu'*it.zu)/(2*n + 2*r);

end