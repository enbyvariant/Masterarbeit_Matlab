function [it, iterations, data] = Interior_Points_LP_alt(iter, nlp,dim)
    
    % Initial values
    n = dim.n;
    m = dim.m;
    r = dim.r;
    en = ones(n,1);
    er = ones(r,1);
    eta = 0.995;

    % Compute starting point
    [it] = Init_LP(nlp,dim);
    iterations = 0;
    data = [];
    p = cutest_setup;
    help = helpers(dim,nlp,it);

    name = genvarname(string(iterations));
    data.(name) = it_log_lp(iterations, it, nlp, 0, 0, 0, help, dim);
    data.success = 0;

    while 1
        if iterations > iter
            break;
        end
        iterations = iterations + 1;

            % update helping variables
            help = helpers(dim,nlp,it);
            
            %set up the affine equation system
            mat = [help.PHI     nlp.A'         nlp.C';
                    nlp.A    sparse(m,m)   sparse(m,r);
                    nlp.C    sparse(r,m)   -help.PSI1   ];
            % disp(rank(full(mat)));
            % disp(size(mat));
            omega_1 = - nlp.c + nlp.A'*it.y + nlp.C'*(it.zl - it.zu) - it.bound_xl*help.Sl1*help.Wl*help.beta_l + it.bound_xu*help.Su1*help.Wu*help.beta_u;
            omega_2 = -(nlp.A*it.x -nlp.b);
            xi = -it.zl - it.bound_cl*help.Tl1*help.Zl*help.rho_l + it.zu + it.bound_cu*help.Tu1*help.Zu*help.rho_u;
            omega_3 = help.PSI1*xi;
            omega = [omega_1; omega_2; omega_3];

            % calculate affine step
            sol = mat\omega;
            if(any(isnan(sol)))
                sol = lsqminnorm(mat,omega);
            end
            del_x_a = sol(1:n);
    
            del_sl_a = it.bound_xl*(del_x_a + help.beta_l);
            del_su_a = it.bound_xu*(-del_x_a + help.beta_u);
            del_wl_a = it.bound_xl*(-it.wl -help.Sl1*help.Wl*del_sl_a);
            del_wu_a = it.bound_xu*(-it.wu - help.Su1*help.Wu*del_su_a);
        
            del_tl_a = it.bound_cl*(nlp.C*del_x_a + help.rho_l);
            del_tu_a = it.bound_cu*(-nlp.C*del_x_a + help.rho_u);
            del_zu_a = it.bound_cu*(-it.zu -help.Tu1*help.Zu*del_tu_a);
            del_zl_a = it.bound_cl*(-it.zl - help.Tl1*help.Zl*del_tl_a);
            
            % calculate duality measure
            mu = (it.sl'*it.wl + it.su'*it.wu + it.tl'*it.zl + it.tu'*it.zu)/(2*(n+r));
            
            % calculate primal affine step length

            step = [del_sl_a; del_su_a; del_tl_a; del_tu_a];
            curr = [it.sl; it.su; it.tl; it.tu];
            [alpha_pri] = step_length(step, curr, 1);

            % calculate dual affine step length
            step = [del_wl_a; del_wu_a; del_zl_a; del_zu_a];
            curr = [it.wl; it.wu; it.zl; it.zu];
            [alpha_dual] = step_length(step, curr, 1);
            
            %------------------------
            % calculate affine duality measure
            sl_a = it.sl + alpha_pri*del_sl_a;
            su_a = it.su + alpha_pri*del_su_a;
            tl_a = it.tl + alpha_pri*del_tl_a;
            tu_a = it.tu + alpha_pri*del_tu_a;
            wl_a = it.wl + alpha_dual*del_wl_a;
            wu_a = it.wu + alpha_dual*del_wu_a;
            zl_a = it.zl + alpha_dual*del_zl_a;
            zu_a = it.zu + alpha_dual*del_zu_a;
            mu_aff = (sl_a'*wl_a + su_a'*wu_a + tl_a'*zl_a + tu_a'*zu_a)/(2*(n+r));
    
            % set the centering parameter sigma
            sigma = (mu_aff/mu)^3;
    
            % set up equation system
            psil = it.bound_cl*(it.zl - sigma*mu*help.Tl1*er + help.Tl1*diag(del_tl_a)*diag(del_zl_a)*er);
            psiu = it.bound_cu*(it.zu - sigma*mu*help.Tu1*er + help.Tu1*diag(del_tu_a)*diag(del_zu_a)*er);
            phil = it.bound_xl*(it.wl - sigma*mu*help.Sl1*en + help.Sl1*diag(del_sl_a)*diag(del_wl_a)*en);
            phiu = it.bound_xu*(it.wu - sigma*mu*help.Su1*en + help.Su1*diag(del_su_a)*diag(del_wu_a)*en);
            
            omega_1 = -nlp.c + nlp.A'*it.y + nlp.C'*(it.zl-it.zu) + it.bound_xl*it.wl - it.bound_xu*it.wu - it.bound_xl*(phil + help.Sl1*help.Wl*help.beta_l) + it.bound_xu*(phiu + help.Su1*help.Wu*help.beta_u);
            xi = -psil + psiu - it.bound_cl*help.Tl1*help.Zl*help.rho_l + it.bound_cu*help.Tu1*help.Zu*help.rho_u;
            omega_3 = help.PSI1*xi;
            omega = [omega_1;omega_2;omega_3];

            % calculate new step
            sol = mat\omega;
            if(any(isnan(sol)))
                sol = lsqminnorm(mat,omega);
            end
            del_x = sol(1:n);
            del_y = -sol(n+1:n+m);
            del_sl = it.bound_xl*(del_x + help.beta_l);
            del_su = it.bound_xu*(-del_x + help.beta_u);
            del_wl = it.bound_xl*(-phil -help.Sl1*help.Wl*del_sl);
            del_wu = it.bound_xu*(-phiu - help.Su1*help.Wu*del_su);
        
            del_tl = it.bound_cl*(nlp.C*del_x + help.rho_l);
            del_tu = it.bound_cu*(-nlp.C*del_x + help.rho_u);
            del_zu = it.bound_cu*(-psiu -help.Tu1*help.Zu*del_tu);
            del_zl = it.bound_cl*(-psil - help.Tl1*help.Zl*del_tl);
            
            %----------------------------
            % calculate primal affine step length

            step = [del_sl; del_su; del_tl; del_tu];
            curr = [it.sl; it.su; it.tl; it.tu];
            [alpha_pri] = step_length(step, curr, eta);

            % calculate dual affine step length
            step = [del_wl; del_wu; del_zl; del_zu];
            curr = [it.wl; it.wu; it.zl; it.zu];
            [alpha_dual] = step_length(step, curr, eta);
            
            % update iterate
            it.x = it.x + alpha_pri*del_x;
            it.sl = it.sl + alpha_pri*del_sl;
            it.su = it.su + alpha_pri*del_su;
            it.tl = it.tl + alpha_pri*del_tl;
            it.tu = it.tu + alpha_pri*del_tu;
            it.y = it.y + alpha_dual*del_y;
            it.zl = it.zl + alpha_dual*del_zl;
            it.zu = it.zu + alpha_dual*del_zu;
            it.wl = it.wl + alpha_dual*del_wl;
            it.wu = it.wu + alpha_dual*del_wu;

            it.mu = (it.sl'*it.wl + it.su'*it.wu + it.tl'*it.zl + it.tu'*it.zu)/(2*(n+r));
            %obj = cutest_obj(it.x);
            name = genvarname(string(iterations));
            data.(name) = it_log_lp(iterations, it, nlp, sigma, alpha_pri, alpha_dual, help, dim);
            convergence = max([data.(name).pri_fea, data.(name).dual_fea, data.(name).compl]);
            if convergence < 10^(-6)
                data.success = 1;
                break
            end
        if it.mu < 1*10^(-12)
            data.success = 1;
            break;
        end
        if it.mu > 1*10^50
            fprintf(1,"did not converge\n");
            break
        end
    end
    
    disp(data);
    cutest_terminate
end