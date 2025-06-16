function [x, sl, su, tl, tu, y, wl, wu, zl, zu, bound_xl,bound_xu, bound_cl, bound_cu] = Interior_Points_Init(H, A, b, c, xl, xu, cl, cu, C)
% Calculate initial iterate of Interior Point QP Method 

    n = size(A, 2);
    m = size(A, 1);
    r = size(C,1);
    en = ones(n,1);

    wl = ones(n,1);
    wu = ones(n,1);

    zl = ones(r,1);
    zu = ones(r,1);

    x = ones(n,1);
    y = rand(m,1);

    % cancel lines where x-component has no upper or no lower bound
    index_l = ones(1,n);
    index_u = ones(1,n);
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

    sl = bound_xl*ones(n,1);
    su = bound_xu*ones(n,1);
    Sl1 = zeros(n);
    Su1 = zeros(n);
    for i = 1:n
        if sl(i) ~= 0
            Sl1(i,i) = 1/sl(i);
        end
        if su(i) ~= 0
            Su1(i,i) = 1/su(i);
        end
    end
    Wl = diag(wl);
    Wu = diag(wu);

    % % Differentiate between the case where inequality constraints exist
    % % or not
     if r == 0
        bound_cl = 0;
        bound_cu = 0;
    
        % helping variables
        betal = bound_xl*(x-xl-sl);
        betau = bound_xu*(-x+xu-su);
        Hphi = H + Sl1*Wl + Su1*Wu;

        % % calculate the affine step with given iterate
        mat = [Hphi A'; A zeros(m)];
        omega_1 = -H*x - c + A'*y - bound_xl*Sl1*Wl*betal + bound_xu*Su1*Wu*betau;
        omega_2 = -(A*x -b);
        omega = [omega_1; omega_2];
        sol = mat\omega;
        del_x = sol(1:n);
        %del_y = sol(n+1:n+m);

        % calculate all variables of expanded system
        del_sl = bound_xl*(del_x + betal);
        del_su = bound_xu*(-del_x + betau);
        del_wl = bound_xl*(-wl - bound_xl*Sl1*Wl*del_sl);
        del_wu = bound_xu*(-wu - bound_xu*Su1*Wu*del_su);

        % ensure positivity s and w of the starting iterate
        sl = bound_xl*max(abs(sl + del_sl),en);
        su = bound_xu*max(abs(su + del_su),en);

        wl = bound_xl*max(abs(wl + del_wl),en);
        wu = bound_xu*max(abs(wu + del_wu),en);

    else
        index_cl = ones(1,r);
        index_cu = ones(1,r);
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
        Zl = diag(zl);
        Zu = diag(zu);
        index_cl = find(index_cl);
        index_cu = find(index_cu);
        bound_cl = sparse(index_cl, index_cl, ones(size(index_cl)), r, r);
        bound_cu = sparse(index_cu, index_cu, ones(size(index_cu)), r, r);

        tl = bound_cl*ones(r,1);
        tu = bound_cu*ones(r,1);
        Tl1 = zeros(r);
        Tu1 = zeros(r);
        for i = 1:r
            if tl(i) ~= 0
            Tl1(i,i) = 1/tl(i);
            end
            if tu(i) ~= 0
                Tu1(i,i) = 1/tu(i);
            end
        end

        % helping variables
        rhol = bound_cl*(C*x-cl-tl);
        rhou = bound_cu*(-C*x + cu - tu);
        betal = bound_xl*(x-xl-sl);
        betau = bound_xu*(-x+xu-su);
        Hphi = H + Sl1*Wl + Su1*Wu;
        psi = Tl1*Zl +Tu1*Zu;
        xi = -bound_cl*(zl + Tl1*Zl*rhol) + bound_cu*(zu + Tu1*Zu*rhou);
        foo = psi\xi;

        % calculate affine step with reduced system
        mat = [Hphi A' C';
            A zeros(m,m + r);
            C zeros(r,m) inv(psi)];
        omega = [-H*x - c + A'*y + C'*bound_cl*zl - C'*bound_cu*zu - bound_xl*Sl1*Wl*betal + bound_xu*Su1*Wu*betau;
            -(A*x-b);
            foo];
        sol = mat\omega;

        del_x = sol(1:n);
        del_z = -sol(n+m+1:n+m+r);

        % calculate all variables of expanded system
        del_sl = bound_xl*(del_x + betal);
        del_su = bound_xu*(-del_x + betau);
        del_wl = bound_xl*(-wl -Sl1*Wl*del_sl);
        del_wu = bound_xu*(-wu - Su1*Wu*del_su);

        del_tl = bound_cl*(C*del_x + rhol);
        del_tu = bound_cu*(-C*del_x + rhou);
        del_zl = -bound_cl*(Tl1*Zl*del_tl + zl);
        del_zu = bound_cu*(-del_z + del_zl);


        % ensure positivity of s, w, t and z (interiority)
        sl = bound_xl*max(abs(sl + del_sl),en);
        su = bound_xu*max(abs(su + del_su),en);

        wl = bound_xl*max(abs(wl + del_wl),en);
        wu = bound_xu*max(abs(wu + del_wu),en);

        tl = bound_cl*max(abs(tl + del_tl),ones(r,1));
        tu = bound_cu*max(abs(tu + del_tu),ones(r,1));

        zl = bound_cl*max(abs(zl + del_zl),ones(r,1));
        zu = bound_cu*max(abs(zu + del_zu),ones(r,1));

    end
    
end