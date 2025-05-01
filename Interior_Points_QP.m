% Implementation of the Interior Points Method for LPs
% jetzt auch in github

function [x,y,zl,zu,wl, wu, sl,su, tl, tu, mu] = Interior_Points_QP(iter, A, b, xl, xu, H, c, cl, cu, C)
    
    % Initial values
    n = size(A, 2);
    m = size(A, 1);
    r = size(C,1);
    en = ones(n,1);
    er = ones(r);
    tao = .5;

    % Compute starting point
    [x, sl, su, tl, tu, y, wl, wu, zl, zu, bound_xl,bound_xu, bound_cl, bound_cu] = Interior_Points_Init(H, A, b, c, xl, xu, cl, cu, C);
    
    for i = 1:iter
        Sl1 = zeros(n);
        Su1 = zeros(n);
        for i = 1:n
            Sl1(i) = max(1/sl(i),0);
            Su1(i) = max(1/su(i),0);
        end
        Wl = diag(wl);
        Wu = diag(wu);

        Tl1 = zeros(r);
        Tu1 = zeros(r);
        for i = 1:r
            Tl1(i) = max(1/tl(i),0);
            Tu1(i) = max(1/tu(i),0);
        end
        Zl = diag(zl);
        Zu = diag(zu);
        
        % calculate affine step
        Hphi = H + Sl1*Wl + Su1*Wu;
        psi = Tl1*Zl + Tu1*Zu;
        rhol = (C*x - cl -tl)*bound_cl;
        rhou = (-C*x + cu -tu)*bound_cu;
        betal = bound_xl*(x -xl -sl);
        betau = bound_xu*(-x + xu -su);
        xi = zl - zu - Tu1*Zu*(rhol+rhou);
        foo = (xi\psi)';
        mat = [Hphi A' C';
            A zeros(m,m+r);
            C zeros(r,m) inv(psi)];
        if r == 0
            bar_l = zeros(size(c));
            bar_u = zeros(size(c));
            omega_3 = zeros(0,1);
        else
            bar_l = C'*bound_cl*zl;
            bar_u = C'*bound_cu*zu;
            omega_3 = rhol + foo;
        end
        omega_1 = H*x + c - A'*y - bar_l + bar_u +Sl1*Wl*betal -Su1*Wu*betau;
        omega_2 = A*x-b;

        omega = [omega_1;omega_2;omega_3];
        sol = (-omega\mat)';

        del_x_a = sol(1:n);
        del_z_a = sol(n+m+1:n+m+r);
        del_sl_a = bound_xl*(del_x_a + betal);
        del_su_a = bound_xu*(-del_x_a + betau);
        del_wl_a = bound_xl*(-wl -Sl1*Wl*del_sl_a);
        del_wu_a = bound_xu*(-wu - Su1*Wu*del_su_a);
    
        del_tl_a = (-del_z_a -xi)\psi;
        del_tu_a = -del_tl_a + rhou + rhol;
        del_zu_a = -zu -Tu1*Zu*del_tu_a;
        del_zl_a = del_z_a + del_zu_a;
        
        % calculate duality measure
        sw = sl'*wl + su'*wu;
        if r == 0
            tz = 0;
        else
            tz = tl'*zl + tu'*zu;
        end
        mu = (sw + tz)/(2*n+2*r);
        
        % calculate affine step length
        step = [del_sl_a; del_su_a; del_wl_a; del_wu_a; del_tl_a; del_tu_a; del_zl_a; del_zu_a];
        curr = [sl; su; wl; wu; tl; tu; zl; zu];
        index = find(step < 0);
        if isempty(index)
            alpha = 1;
        else
            alpha = min(curr(index)./(-step(index)));
            if alpha > 1
                alpha = 1;
            end
        end
        
        % calculate affine duality measure
        slwl_a = (sl+ alpha*del_sl_a)' * (wl + alpha*del_wl_a);
        suwu_a = (su + alpha*del_su_a)' * (wu + del_wu_a);
        if r == 0
            tlzl_a = 0;
            tuzu_a = 0;
        else
            tlzl_a = (tl + alpha*del_tl_a)'*(zl + alpha*del_zl_a);
            tuzu_a = (tu + alpha*del_tu_a)'*(zu + alpha*del_zu_a);
        end
        mu_aff = (slwl_a + suwu_a + tlzl_a + tuzu_a)/(2*n+2*r);

        % set the centering parameter sigma
        sigma = (mu_aff/mu)^3;

        % calculate new step
        xi = zl - sigma*mu*Tl1*er + Tl1*diag(del_tl_a)*diag(del_zl_a)*er - zu + sigma*mu*Tu1*er - Tu1*diag(del_tu_a)*diag(del_zu_a)*er - Tu1*Zu*(rhol+rhou);
        foo = (xi\psi)';
        omega_1 = H*x + c - A'*y - bar_l + bar_u +Sl1*Wl*betal -Su1*Wu*betau - sigma*mu*(Sl1*en) + sigma*mu*(Su1*en) + Sl1*diag(del_sl_a)*diag(del_wl_a)*en - Su1*diag(del_su_a)*diag(del_wu_a)*en;
        if r ~= 0
            omega_3 = rhol + foo;
        end
        omega = [omega_1;omega_2;omega_3];
        sol = (-omega\mat)';
        del_x = sol(1:n);
        del_y = sol(n+1:n+m);
        del_z = sol(n+m+1:n+m+r);
        del_sl = bound_xl*(del_x + betal);
        del_su = bound_xu*(-del_x + betau);
        del_wl = bound_xl*(-wl + sigma*mu*Sl1*en - Sl1*diag(del_sl_a)*diag(del_wl_a)*en -Sl1*Wl*del_sl);
        del_wu = bound_xu*(-wu + sigma*mu*Su1*en - Su1*diag(del_su_a)*diag(del_wu_a)*en - Su1*Wu*del_su);
    
        del_tl = (-del_z -xi)\psi;
        del_tu = -del_tl + rhou + rhol;
        del_zu = -zu + sigma*mu*Tu1*er - Tu1*diag(del_tu_a)*diag(del_zu_a)*er -Tu1*Zu*del_tu;
        del_zl = del_z + del_zu;
        
        % calculate step length alpha
        step = [del_sl; del_su; del_wl; del_wu; del_tl; del_tu; del_zl; del_zu];
        curr = [sl; su; wl; wu; tl; tu; zl; zu];
        index = find(step < 0);
        if isempty(index)
            alpha = 1;
        else
            alpha = min(tao*curr(index)./(-step(index)));
            if alpha > 1
                alpha = 1;
            end
        end
        
        % update iterate
        x = x + alpha*del_x;
        y = y + alpha*del_y;
        zl = zl + alpha*del_zl;
        zu = zu + alpha*del_zu;
        wl = wl + alpha*del_wl;
        wu = wu + alpha*del_wu;
        sl = sl + alpha*del_sl;
        su = su + alpha*del_su;
        tl = tl + alpha*del_tl;
        tu = tu + alpha*del_tu;

        tao = tao/2;
    end

end