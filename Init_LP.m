function [x, sl, su, tl, tu, y, wl, wu, zl, zu, bound_xl, bound_xu, bound_cl, bound_cu] = Init_LP(A, C, cl, cu, xl, xu, b, c, dim)
%
% 
    % Helping vectors
    en = ones(dim.n,1);
    er = ones(dim.r,1);
    
    index_l = ones(1,dim.n);
    index_u = ones(1,dim.n);
    wl = ones(dim.n, 1);
    wu = ones(dim.n, 1);
    for i = 1:dim.n
        if xl(i) < -10^3
            xl(i) = 0;
            wl(i) = 0; 
            index_l(i) = 0;
        end
        if xu(i) > 10^3
            xu(i) = 0;
            wu(i) = 0;
            index_u(i) = 0;
        end
    end
    index_l = find(index_l);
    index_u = find(index_u);
    bound_xl = sparse(index_l, index_l, ones(size(index_l)), dim.n, dim.n);
    bound_xu = sparse(index_u, index_u, ones(size(index_u)), dim.n, dim.n);

    %bound_xu = zeros(n);
    
    index_cl = ones(1, dim.r);
    index_cu = ones(1, dim.r);
    zl = ones(dim.r, 1);
    zu = ones(dim.r, 1);
    for i = 1:dim.r
        if cl(i) < -10^3
            cl(i) = 0;
            zl(i) = 0;
            index_cl(i) = 0;
        end
        if cu(i) > 10^3
            cu(i) = 0;
            zu(i) = 0;
            index_cu(i) = 0;
        end
    end
    index_cl = find(index_cl);
    index_cu = find(index_cu);
    bound_cl = sparse(index_cl, index_cl, ones(size(index_cl)), dim.r, dim.r);
    bound_cu = sparse(index_cu, index_cu, ones(size(index_cu)), dim.r, dim.r);

    % Compute starting point
    M = [A zeros(dim.m, dim.n + dim.r);
        2*eye(dim.n) -eye(dim.n) zeros(dim.n, dim.r);
        2*C zeros(dim.r, dim.n) -eye(dim.r)];
    vector = [b; bound_xl*xl + bound_xu*xu; bound_cl*cl + bound_cu*cu]; 
    temp = (M*M'\vector);
    a_prim = M' * temp;
    x = a_prim(1:dim.n);
    sl = bound_xl*(x - xl);
    su = bound_xu*(-x + xu);
    tl = bound_cl*(C*x - cl);
    tu = bound_cu*(-C*x + cu);
    
    N = [A' -eye(dim.n) C'];
    a_dual = N' * (N*N'\c);
    y = a_dual(1:dim.m);
    zl = a_dual(dim.m + dim.n + 1 : dim.m + dim.n + dim.r);
    zu = zeros(dim.r, 1);
    wl = a_dual(dim.m + 1 : dim.m + dim.n);
    wu = zeros(dim.n, 1);
    
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
        x = ones(dim.n,1);
    end
    if any(isnan(sl))
        sl = bound_xl*ones(dim.n,1);
    end
    if any(isnan(su))    
        su = bound_xu*ones(dim.n,1);
    end
    if any(isnan(tl))

        tl = bound_cl*ones(dim.r,1);
    end
    if any(isnan(tu))
        tu = bound_cu*ones(dim.r,1);
    end
    if any(isnan(wl))
        wl = bound_xl*ones(dim.n,1);
    end

    if any(isnan(wu))
        wu = bound_xu*ones(dim.n,1);
    end
    if any(isnan(zl))
        zl = bound_cl*ones(dim.r,1);
    end
    if any(isnan(zu))
        zu = bound_cu*ones(dim.r,1);
    end
    if any(isnan(y))
        y = ones(dim.m,1);
    end

end