function [x, sl, su, tl, tu, y, wl, wu, zl, zu, bound_xl, bound_xu, bound_cl, bound_cu] = Init_LP(A, C, cl, cu, xl, xu, b, c)
%
% 
    % Initial values
    n = size(A, 2);
    m = size(A, 1);
    r = size(C, 1);
    en = ones(n,1);
    er = ones(r,1);
    
    index_l = ones(1,n);
    index_u = ones(1,n);
    wl = ones(n,1);
    wu = ones(n,1);
    for i = 1:n
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
    bound_xl = sparse(index_l, index_l, ones(size(index_l)), n, n);
    bound_xu = sparse(index_u, index_u, ones(size(index_u)), n, n);

    %bound_xu = zeros(n);
    
    index_cl = ones(1,r);
    index_cu = ones(1,r);
    zl = ones(r,1);
    zu = ones(r,1);
    for i = 1:r
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
    bound_cl = sparse(index_cl, index_cl, ones(size(index_cl)), r, r);
    bound_cu = sparse(index_cu, index_cu, ones(size(index_cu)), r, r);

    % Compute starting point
    M = [A zeros(m,n+r);2*bound_xl*eye(n) -eye(n) zeros(n,r);2*bound_cl*C zeros(r,n) -eye(r)];
    vector = [b; bound_xl*xl + bound_xu*xu; bound_cl*cl + bound_cu*cu];

    a_prim = M' * (M*M'\vector);
    x = a_prim(1:n);
    sl = bound_xl*(x - xl);
    su = bound_xu*(-x + xu);
    tl = bound_cl*(C*x - cl);
    tu = bound_cu*(-C*x + cu);
    
    N = [A' eye(n) C'];
    a_dual = N' * (N*N'\c);
    y = a_dual(1:m);
    zl = a_dual(m+n+1:m+n+r);
    zu = zeros(r,1);
    wl = a_dual(m+1:m+n);
    wu = zeros(n,1);
    
    delta_pri = 3/2*max([-sl;-su;-wl;-wu;0]);
    delta_dual = -3/2*max([wl;wu;zl;zu; 0]);
    
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
end