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
        it.mu = (it.sl'*it.wl + it.su'*it.wu + it.tl'*it.zl + it.tu'*it.zu)/(2*m + 2*r);
        if it.mu < 10^(-15)
            break
        end

       % calculate step direction
        [nlp] = cutest_iterate(it, nlp, dim,p);
        [help] = helpers(dim, nlp, it);

        phi_l = it.wl - help.Sl1*it.mu*en;
        phi_u = it.wu - help.Su1*it.mu*en;
        psi_l = it.zl - help.Tl1*it.mu*er;
        psi_u = it.zu - help.Tu1*it.mu*er;

        mat = [nlp.H + help.PHI   nlp.A'  nlp.C';
             nlp.A   zeros(m)    zeros(m,r);
             nlp.C  zeros(r,m)   -help.PSI1];
        nu = - nlp.grad - phi_l -help.Sl1*help.Wl*help.beta_l ...
            + phi_u + help.Su1*help.Wu*help.beta_u;
        xi = -psi_l - help.Tl1*help.Zl*help.rho_l ...
        + psi_u + help.Tu1*help.Zu*help.rho_u;
        omega = [nu; -nlp.c_e ; help.PSI1*xi]
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