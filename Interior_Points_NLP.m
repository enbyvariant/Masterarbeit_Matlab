function [nlp,it, iter] = Interior_Points_NLP(max_iter)
% Interior Point Method for general NLPs
    
    % Initialization
    [it,nlp,dim] = Interior_gen_Init();
    n = dim.n;
    m = dim.m;
    r = dim.r;
    en = ones(n,1);
    er = ones(r,1);

    eta = .995;

    p = cutest_setup();

    iter = 0;
    
    while 1
        if iter > max_iter
            break
        else
            iter = iter + 1;
        end


        % duality measure
        it.mu = (it.sl'*it.wl + it.su'*it.wu + it.tl'*it.zl + it.tu'*it.zu)/(2*n + 2*r);
        if it.mu < 10^(-15)
            break
        end
    
        % calculate affine step
        [nlp] = cutest_iterate(it, nlp, dim,p);
        [help] = helpers_nlp(dim, nlp, it,p);


        %disp(help.PHI)

        mat = [nlp.H + help.PHI   nlp.A'  nlp.C';
             nlp.A   zeros(m)    zeros(m,r);
             nlp.C  zeros(r,m)   -help.PSI1];

        nu = - nlp.grad + nlp.A'*it.y + nlp.C'*(it.zl - it.zu) -it.bound_xl*help.Sl1*help.Wl*help.beta_l ...
            + it.bound_xu*help.Su1*help.Wu*help.beta_u;
        xi = -it.zl - it.bound_cl*help.Tl1*help.Zl*help.rho_l ...
        + it.zu + it.bound_cu*help.Tu1*help.Zu*help.rho_u;
        omega = [nu; -nlp.c_e; help.PSI1*xi];

        sol = mat\omega;
        del_x = sol(1:n);

        del_sl = it.bound_xl*(del_x + help.beta_l);
        del_su = it.bound_xu*(-del_x + help.beta_u);
        del_wl = it.bound_xl*(-it.wl -help.Sl1*help.Wl*del_sl);
        del_wu = it.bound_xu*(-it.wu - help.Su1*help.Wu*del_su);
        
        del_tl = it.bound_cl*(nlp.C*del_x + help.rho_l);
        del_tu = it.bound_cu*(-nlp.C*del_x + help.rho_u);
        del_zl = it.bound_cl*(-it.zl -help.Tl1*help.Zl*del_tl);
        del_zu = it.bound_cu*(-it.zu - help.Tu1*help.Zu*del_tu);

        % calculate step length
        step = [del_sl; del_su; del_wl; del_wu; del_tl; del_tu; del_zl; del_zu];
        curr = [it.sl; it.su; it.wl; it.wu; it.tl; it.tu; it.zl; it.zu];
        [alpha] = step_length(step, curr, 1);

        % compute affine duality measure
        sl_a = it.sl + alpha*del_sl;
        su_a = it.su + alpha*del_su;
        tl_a = it.tl + alpha*del_tl;
        tu_a = it.tu + alpha*del_tu;
        wl_a = it.wl + alpha*del_wl;
        wu_a = it.wu + alpha*del_wu;
        zl_a = it.zl + alpha*del_zl;
        zu_a = it.zu + alpha*del_zu;

        mu_aff = (sl_a'*wl_a + su_a'*wu_a + tl_a'*zl_a + tu_a'*zu_a)/(2*(n+r));

        % centering parameter
        sigma = (mu_aff/it.mu)^3;
        
        % calculate step direction
        phi_l = it.bound_xl*(it.wl - help.Sl1*sigma*it.mu*en + help.Sl1*diag(del_sl)*diag(del_wl)*en);
        phi_u = it.bound_xu*(it.wu - help.Su1*sigma*it.mu*en+ help.Su1*diag(del_su)*diag(del_wu)*en);
        psi_l = it.bound_cl*(it.zl - sigma*it.mu*help.Tl1*er+ help.Tl1*diag(del_tl)*diag(del_zl)*er);
        psi_u = it.bound_cu*(it.zu - sigma*it.mu*help.Tu1*er+ help.Tu1*diag(del_tu)*diag(del_zu)*er);
        
        nu = - nlp.grad + nlp.A'*it.y + nlp.C'*(it.zl-it.zu) + it.bound_xl*it.wl - it.bound_xu*it.wu ...
                - it.bound_xl*(phi_l + help.Sl1*help.Wl*help.beta_l) + it.bound_xu*(phi_u + help.Su1*help.Wu*help.beta_u);
        xi = -psi_l - it.bound_cl*help.Tl1*help.Zl*help.rho_l ...
        + psi_u + it.bound_cu*help.Tu1*help.Zu*help.rho_u;
        omega = [nu; -nlp.c_e ; help.PSI1*xi];

        disp(it.mu);
    
        sol = mat\omega;
        del_x = sol(1:n);
        del_y = - sol(n+1:n+m);

        del_sl = it.bound_xl*(del_x + help.beta_l);
        del_su = it.bound_xu*(-del_x + help.beta_u);
        del_wl = it.bound_xl*(-phi_l -help.Sl1*help.Wl*del_sl);
        del_wu = it.bound_xu*(-phi_u - help.Su1*help.Wu*del_su);
        
        del_tl = it.bound_cl*(nlp.C*del_x + help.rho_l);
        del_tu = it.bound_cu*(-nlp.C*del_x + help.rho_u);
        del_zu = it.bound_cu*(-psi_u -help.Tu1*help.Zu*del_tu);
        del_zl = it.bound_cl*(-psi_l - help.Tl1*help.Zl*del_tl);

        % calculate step length alpha
        step = [del_sl; del_su; del_wl; del_wu; del_tl; del_tu; del_zl; del_zu];
        curr = [it.sl; it.su; it.wl; it.wu; it.tl; it.tu; it.zl; it.zu];
        index = find(step < 0);
        if isempty(index)
            alpha = 1;
        else
            alpha = min(eta*curr(index)./(-step(index)));
            if alpha > 1
                alpha = 1;
            end
        end

        % compute new iterate
        it.x = it.x + alpha*del_x;
        %disp(it.x)
        it.y = it.y + alpha*del_y;
        it.sl = it.sl + alpha*del_sl;
        it.su = it.su +alpha*del_su;
        it.wl = it.wl + alpha*del_wl;
        it.wu = it.wu + alpha*del_wu;
        it.tl = it.tl + alpha*del_tl;
        it.tu = it.tu + alpha*del_tu;
        it.zl = it.zl + alpha*del_zl;
        it.zu = it.zu + alpha*del_zu;

    end
    
    cutest_terminate;

end