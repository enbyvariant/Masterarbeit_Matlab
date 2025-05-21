function [x, sl, su, tl, tu, y, wl, wu, zl, zu, bound_xl,bound_xu, bound_cl, bound_cu] = Interior_Points_Init(H, A, b, c, xl, xu, cl, cu, C)
% Calculate initial iterate of Interior Point QP Method 

    n = size(A, 2);
    m = size(A, 1);
    r = size(C,1);
    en = ones(n,1);

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

    % % Differentiate between the case where inequality constraints exist
    % % or not
     if r == 0
        bound_cl = NaN;
        bound_cu = NaN;
    
        %helping variables
        Sl = diag(sl);
        Su = diag(su);
        % betal = bound_xl*(x-xl-sl);
        % betau = bound_xu*(-x+xu-su);
        % Hphi = H + Sl1*Wl + Su1*Wu;
            % calculate affine step
            mat = [H     zeros(n)    zeros(n)      -A'       -eye(n)     eye(n);
                zeros(n)    Wl       zeros(n)   zeros(n,m)   Sl        zeros(n);
                zeros(n) zeros(n)       Wu      zeros(n,m)  zeros(n)    Su;
                   A     zeros(m,n)  zeros(m,n)  zeros(m)  zeros(m,n)  zeros(m,n);
                 eye(n)   -eye(n)    zeros(n)   zeros(n,m)  zeros(n)   zeros(n);
                -eye(n)  zeros(n)    -eye(n)    zeros(n,m)  zeros(n)   zeros(n)];
            om1 = -(H*x + c - A' * y - wl + wu);
            om2 = -Wl*Sl*en;
            om3 = -Wu*Su*en;
            om4 = -(A*x - b);
            om5 = -(x - xl - sl);
            om6 = -(-x + xu -su);
            omega = [om1;om2;om3;om4;om5;om6];
            sol = mat\omega;
            %del_x_a = sol(1:n);
            del_sl = sol(n+1:2*n);
            del_su = sol(2*n+1:3*n);
            del_wl = sol(3*n+m+1:4*n+m);
            del_wu = sol(4*n+m+1:5*n+m);

        % % calculate the affine step with given iterate
        % mat = [Hphi A'; A zeros(m)];
        % omega_1 = -H*x - c + A'*y - bound_xl*Sl1*Wl*betal + bound_xu*Su1*Wu*betau;
        % omega_2 = -(A*x -b);
        % omega = [omega_1; omega_2];
        % sol = mat\omega;
        % del_x = sol(1:n);
        % 
        % % calculate all variables of expanded system
        % del_sl = bound_xl*(del_x + betal);
        % del_su = bound_xu*(-del_x + betau);
        % del_wl = bound_xl*(-wl -Sl1*Wl*del_sl);
        % del_wu = bound_xu*(-wu - Su1*Wu*del_su);

        % ensure positivity s and w of the starting iterate
        sl = bound_xl*max(abs(sl + del_sl),en);
        su = bound_xu*max(abs(su + del_su),en);

        wl = bound_xl*max(abs(wl + del_wl),en);
        wu = bound_xu*max(abs(wu + del_wu),en);
    else
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
        xi = -bound_cl*(zl + Tl1*Zl*rhol) + bound_cu*(zu + Tu1*Zu*rhou);
        foo = xi/psi;

        % calculate affine step with reduced system
        mat = [Hphi A' C';
            A zeros(m,m + r);
            C zeros(r,m) inv(psi)];
        omega = [-H*x - c + A'*y + bound_cl*C'*zl - bound_cu*C'*zu - bound_xl*Sl1*Wl*betal + bound_xu*Su1*Wu*betau;
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