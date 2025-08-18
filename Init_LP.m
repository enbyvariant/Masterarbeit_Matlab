function [it, step_a] = Init_LP(nlp, dim)
%
% 
    % Helping vectors
    en = ones(dim.n,1);
    er = ones(dim.r,1);
    fn = zeros(dim.n,1);
    fr = zeros(dim.r,1);

    % set up bound_xl and bound_xu to discern which bounds exist
    index_l = ones(1,dim.n);
    index_u = ones(1,dim.n);
    it.wl = ones(dim.n, 1);
    it.wu = ones(dim.n, 1);

    nlp.H = zeros(dim.n);

    for i = 1:dim.n
        if nlp.xl(i) < -10^3
            nlp.xl(i) = 0;
            it.wl(i) = 0; 
            index_l(i) = 0;
        end
        if nlp.xu(i) > 10^3
            nlp.xu(i) = 0;
            it.wu(i) = 0;
            index_u(i) = 0;
        end
    end

    index_l = find(index_l);
    index_u = find(index_u);
    it.bound_xl = sparse(index_l, index_l, ones(size(index_l)), dim.n, dim.n);
    it.bound_xu = sparse(index_u, index_u, ones(size(index_u)), dim.n, dim.n);

    % set up bound_xl and bound_xu to discern which bounds exist
    index_cl = ones(1, dim.r);
    index_cu = ones(1, dim.r);
    it.zl = ones(dim.r, 1);
    it.zu = ones(dim.r, 1);

    for i = 1:dim.r
        if nlp.cl(i) < -10^3
            nlp.cl(i) = 0;
            it.zl(i) = 0;
            index_cl(i) = 0;
        end
        if nlp.cu(i) > 10^3
            nlp.cu(i) = 0;
            it.zu(i) = 0;
            index_cu(i) = 0;
        end
    end

    index_cl = find(index_cl);
    index_cu = find(index_cu);
    it.bound_cl = sparse(index_cl, index_cl, ones(size(index_cl)), dim.r, dim.r);
    it.bound_cu = sparse(index_cu, index_cu, ones(size(index_cu)), dim.r, dim.r);

    % Compute primal starting point from constraints
    M = [nlp.A zeros(dim.m, dim.n + dim.r);
        2*eye(dim.n) -eye(dim.n) zeros(dim.n, dim.r);
        2*nlp.C zeros(dim.r, dim.n) -eye(dim.r)];
    vector = [nlp.b; it.bound_xl*nlp.xl + it.bound_xu*nlp.xu; it.bound_cl*nlp.cl + it.bound_cu*nlp.cu];

    temp = (M*M'\vector);
    a_prim = M' * temp;
    it.x = a_prim(1:dim.n);
    %it.x = [1;1];
    it.sl = it.bound_xl*(it.x - nlp.xl);
    it.su = it.bound_xu*(-it.x + nlp.xu);
    it.tl = it.bound_cl*(nlp.C*it.x - nlp.cl);
    it.tu = it.bound_cu*(-nlp.C*it.x + nlp.cu);
    
    % Compute dual starting point from Lagrange condition
    N = [nlp.A' -eye(dim.n) nlp.C'];
    a_dual = N' * (N*N'\nlp.c);
    it.y = a_dual(1:dim.m);
    it.zl = a_dual(dim.m + dim.n + 1 : dim.m + dim.n + dim.r);
    it.zu = zeros(dim.r, 1);
    it.wl = a_dual(dim.m + 1 : dim.m + dim.n);
    it.wu = zeros(dim.n, 1);
    
    % ensure positivity of starting point
    delta_pri = 3/2*max([-it.sl;-it.su;-it.tl;-it.tu;0]);
    delta_dual = 3/2*max([-it.wl;-it.wu;-it.zl;-it.zu; 0]);
    
    it.sl = it.bound_xl*(it.sl + delta_pri * en);
    it.su = it.bound_xu*(it.su + delta_pri * en);
    it.tl = it.bound_cl*(it.tl + delta_pri * er);
    it.tu = it.bound_cu*(it.tu + delta_pri * er);
    it.wl = it.bound_xl*(it.wl + delta_dual * en);
    it.wu = it.bound_xu*(it.wu + delta_dual * en);
    it.zl = it.bound_cl*(it.zl + delta_dual * er);
    it.zu = it.bound_cu*(it.zu + delta_dual * er);
    
    % improve centering of starting point
    num = it.sl' * it.wl + it.su'*it.wu + it.tl'*it.zl + it.tu'*it.zu;
    delta_pri = 1/2*num/(en' * it.wl + en'*it.wu + er'*it.zl + er'*it.zu);
    delta_dual = 1/2*num/(en' * it.sl + en'*it.su + er'*it.tl + er'*it.tu);

    it.sl = it.bound_xl*(it.sl + delta_pri * en);
    it.su = it.bound_xu*(it.su + delta_pri * en);
    it.tl = it.bound_cl*(it.tl + delta_pri * er);
    it.tu = it.bound_cu*(it.tu + delta_pri * er);
    it.wl = it.bound_xl*(it.wl + delta_dual * en);
    it.wu = it.bound_xu*(it.wu + delta_dual * en);
    it.zl = it.bound_cl*(it.zl + delta_dual * er);
    it.zu = it.bound_cu*(it.zu + delta_dual * er);
    
    it.x = ones(dim.n,1);
    it.sl = it.bound_xl*ones(dim.n,1);
    it.su = it.bound_xu*ones(dim.n,1);
    it.tl = it.bound_cl*ones(dim.r,1);
    it.tu = it.bound_cu*ones(dim.r,1);
    it.wl = it.bound_xl*ones(dim.n,1);
    it.wu = it.bound_xu*ones(dim.n,1);
    it.zl = it.bound_cl*ones(dim.r,1);
    it.zu = it.bound_cu*ones(dim.r,1);

    % calculate initial duality measure
    it.mu = (it.sl'*it.wl + it.su'*it.wu + it.tl'*it.zl + it.tu'*it.zu)/(2*(dim.n+dim.r));
    
    % helping variable
    step_a.del_sl = fn;
    step_a.del_su = fn;
    step_a.del_wl = fn;
    step_a.del_wu = fn;
    step_a.del_tl = fr;
    step_a.del_tu = fr;
    step_a.del_zl = fr;
    step_a.del_zu = fr;

    % Failsave, if some linear equation system is impossible to solve
    if any(isnan(it.x))
        it.x = ones(dim.n,1);
    end
    if any(isnan(it.sl))
        it.sl = it.bound_xl*ones(dim.n,1);
    end
    if any(isnan(it.su))    
        it.su = it.bound_xu*ones(dim.n,1);
    end
    if any(isnan(it.tl))
        it.tl = it.bound_cl*ones(dim.r,1);
    end
    if any(isnan(it.tu))
        it.tu = it.bound_cu*ones(dim.r,1);
    end
    if any(isnan(it.wl))
        it.wl = it.bound_xl*ones(dim.n,1);
    end

    if any(isnan(it.wu))
        it.wu = it.bound_xu*ones(dim.n,1);
    end
    if any(isnan(it.zl))
        it.zl = it.bound_cl*ones(dim.r,1);
    end
    if any(isnan(it.zu))
        it.zu = it.bound_cu*ones(dim.r,1);
    end
    if any(isnan(it.y))
        it.y = ones(dim.m,1);
    end

end