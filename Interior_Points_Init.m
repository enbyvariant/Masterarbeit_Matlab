function [x, sl, su, tl, tu, y, wl, wu, zl, zu, bound_xl,bound_xu, bound_cl, bound_cu] = Interior_Points_Init(H, A, b, c, xl, xu, cl, cu, C)
% Calculate initial iterate of Interior Point QP Method 

    n = size(A, 2);
    m = size(A, 1);
    r = size(C,1);
    e = ones(n,1);

    sl = rand(n,1);
    su = rand(n,1);
    wl = rand(n,1);
    wu = rand(n,1);

    tl = rand(r,1);
    tu = rand(r,1);
    zl = rand(r,1);
    zu = rand(r,1);

    x = rand(n,1);
    y = rand(m,1);

    Sl1 = diag(ones(n,1)./sl);
    Wl = diag(wl);
    Su1 = diag(ones(n,1)./su);
    Wu = diag(wu);
    Tl1 = diag(ones(r,1)./tl);
    Zl = diag(zl);
    Tu1 = diag(ones(r,1)./tu);
    Zu = diag(zu);


    % cancel lines where x-component has no upper or no lower bound
    bound_xl = eye(n);
    bound_xu = eye(n);
    for i = n:1
        if xl(i) < -10^3
            xl(i) = 0;
            sl(i) = 0; 
            wl(i) = 0; 
            bound_xl(i,i) = 0;
        end
        if xu > 10^3
            xu(i) = 0;
            su(i) = 0;
            wu(i) = 0;
            bound_xu(i,i) = 0;
        end
    end
    Sl1 = bound_xl*Sl1;
    Wl = bound_xl*Wl;
    Su1 = bound_xu*Su1;
    Wu = bound_xu*Wu;

    %cancel lines where Cx has no upper or no lower bound
    bound_cl = eye(r);
    bound_cu = eye(r);
    for i = 1:r
        if cl(i) < -10^3
            cl(i) = 0;
            tl(i) = 0;
            zl(i) = 0;
            bound_cl(i,i) = 0;
            
        end
        if cu(i) > 10^3
            cu(i) = 0;
            tu(i) = 0;
            zu(i) = 0;
            bound_cu(i,i) = 0;
        end
    end
    Tl1 = bound_cl*Tl1;
    Zl = bound_cl*Zl;
    Tu1 = bound_cu*Tu1;
    Zu = bound_cu*Zu;

    
    % helping variables
    rhol = bound_cl*(C*x-cl-tl);
    rhou = bound_cu*(-C*x + cu - tu);
    betal = bound_xl*(x-xl-sl);
    betau = bound_xu*(-x+xu-su);
    Hphi = H + Sl1*Wl + Su1*Wu;
    psi = Tl1*Zl +Tu1*Zu;
    xi = zl-zu-Tu1*Zu*(rhol + rhou);
    foo = xi/psi;

    % calculate affine step with reduced system
    mat = [Hphi A' C';
        A zeros(m,m + r);
        C zeros(r,m) inv(psi)];
    omega = [H*x + c - A'*y - C'*zl + C'*zu;
        A*x-b;
        rhol + foo];
    sol = (-omega\mat)';
    
    del_x = sol(1:n);
    del_z = sol(n+m+1:n+m+r);

    % calculate all variables of expanded system
    del_sl = bound_xl*(del_x + betal);
    del_su = bound_xu*(-del_x + betau);
    del_wl = bound_xl*(-wl -Sl1*Wl*del_sl);
    del_wu = bound_xu*(-wu - Su1*Wu*del_su);

    del_tl = (-del_z -xi)\psi;
    del_tu = -del_tl + rhou + rhol;
    del_zu = -zu -Tu1*Zu*del_tu;
    del_zl = del_z + del_zu;

    
    % ensure positivity of s, w, t and z (interiority)
    sl = bound_xl*max(abs(sl + del_sl),e);
    su = bound_xu*max(abs(su + del_su),e);

    wl = bound_xl*max(abs(wl + del_wl),e);
    wu = bound_xu*max(abs(wu + del_wu),e);

    tl = bound_cl*max(abs(tl + del_tl),ones(r,1));
    tu = bound_cu*max(abs(tu + del_tu),ones(r,1));

    zl = bound_cl*max(abs(zl + del_zl),ones(r,1));
    zu = bound_cu*max(abs(zu + del_zu),ones(r,1));
end