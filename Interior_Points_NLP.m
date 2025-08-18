function [x,y,sl,su,wl,wu,tl, tu,zl,zu,obj] = Interior_Points_NLP(max_iter)
% Interior Point Method for general NLPs

    
    % Initialization
    [it, dim, nlp] = Interior_gen_Init();
    en = ones(n,1);
    er = ones(r,1);

    iter = 0;
    
    while 1
        if iter > max_iter
            break
        else
            iter = iter + 1;
        end

        % duality measure
        it.mu = (sl'*wl + su'*wu + tl'*zl + tu'*zu)/(2*m + 2*r);
        if it.mu < 10^(-15)
            break
        end
        % calculate variables for new iterate
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
        Zl = diag(zl);
        Zu = diag(zu);

        PHI = Sl1*Wl+Su1*Wu;
        PSI = Tl1*Zl+Tu1*Zu;
        PSI1 = zeros(r);
        for i = 1:r
            PSI1(i,i) = 1/PSI(i,i);
        end
        [help] = helpers(dim, nlp, it);
        % calculate step direction
        [grad, H, cons_e, cons_i, cons_e_grad, C, obj] = cutest_iterate(x,y,zl,zu,equ,inequ);
        beta_l = x - xl - sl;
        beta_u = -x + xu - su;
        rho_l = cons_i - tl;
        rho_u = -cons_i - tu;
        phi_l = wl - Sl1*mu*en;
        phi_u = wu - Su1*mu*en;
        psi_l = zl - Tl1*mu*er;
        psi_u = zu - Tu1*mu*er;

        mat = [help.Hphi   cons_e_grad'  C';
             cons_e_grad   zeros(m)    zeros(m,r);
             C  zeros(r,m)   -PSI1];
        nu = - grad - phi_l -Sl1*Wl*beta_l + phi_u + Su1*Wu*beta_u;
        xi = -psi_l - Tl1*Zl*rho_l + psi_u + Tu1*Zu*rho_u;
        omega = [nu; -cons_e ; PSI1*xi];
        sol = mat\omega;
        del_x = sol(1:n);
        del_y = - sol(n+1:n+m);
        del_z = - sol(n+m+1:n+m+r);
        del_sl = bound_xl*(del_x + beta_l);
        del_su = bound_xu*(-del_x + beta_u);
        del_wl = bound_xl*(-phi_l -Sl1*Wl*del_sl);
        del_wu = bound_xu*(-phi_u - Su1*Wu*del_su);
        
        del_tl = bound_cl*(C*del_x + rho_l);
        del_tu = bound_cu*(-C*del_x + rho_u);
        del_zu = bound_cu*(-psi_u -Tu1*Zu*del_tu);
        del_zl = bound_cl*(-psi_l - Tl1*Zl*del_tl);

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

        % compute new iterate
        x = x + alpha*del_x;
        y = y + alpha*del_y;
        sl = sl + alpha*del_sl;
        su = su +alpha*del_su;
        wl = wl + alpha*del_wl;
        wu = wu + alpha*del_wu;
        tl = tl + alpha*del_tl;
        tu = tu + alpha*del_tu;
        zl = zl + alpha*del_zl;
        zu = zu + alpha*del_zu;


    end

end