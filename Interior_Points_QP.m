% Implementation of the Interior Points Method for LPs
% jetzt auch in github

function [it, mu_n, obj, iterations, data] = Interior_Points_QP(iter, nlp, dim)
    
    % Initial values
    n = dim.n;
    m = dim.m;
    r = dim.r;
    en = ones(n,1);
    er = ones(r,1);
    eta = 0.995;

    % Compute starting point
    [it] = Interior_Points_Init(nlp,dim);
    %[x, sl, su, tl, tu, y, wl, wu, zl, zu, bound_xl,bound_xu, bound_cl, bound_cu] = Interior_Points_Init(H, A, b, c, xl, xu, cl, cu, C,nlp,dim);
    iterations = 0;

    data = [];
    p = cutest_setup;


    while 1
        if iterations > iter
            break;
        end
        iterations = iterations + 1;


            % update helping variables
            [help] = helpers(dim, nlp, it);
            
            % calculate affine step
            Hphi = nlp.H + help.PHI;

             mat = [Hphi     nlp.A'         nlp.C';
                   nlp.A   zeros(m)   zeros(m,r);
                   nlp.C   zeros(r,m)   -help.PSI1 ];
            omega_1 = -nlp.H*it.x - nlp.c + nlp.A'*it.y + nlp.C'*(it.zl - it.zu) - ...
                it.bound_xl*help.Sl1*help.Wl*help.beta_l + it.bound_xu*help.Su1*help.Wu*help.beta_u;
            omega_2 = -(nlp.A*it.x -nlp.b);
            xi = -it.zl - it.bound_cl*help.Tl1*help.Zl*help.rho_l + it.zu + it.bound_cu*help.Tu1*help.Zu*help.rho_u;
            omega_3 = help.PSI1*xi;
            omega = [omega_1; omega_2; omega_3];
            sol = mat\omega;
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
            
            % calculate affine step length
            step = [del_sl_a; del_su_a; del_wl_a; del_wu_a;...
                del_tl_a; del_tu_a; del_zl_a; del_zu_a];
            curr = [it.sl; it.su; it.wl; it.wu; it.tl; it.tu; it.zl; it.zu];
            [alpha] = step_length(step, curr, 1);
            
            % calculate affine duality measure
            sl_a = it.sl + alpha*del_sl_a;
            su_a = it.su + alpha*del_su_a;
            wl_a = it.wl + alpha*del_wl_a;
            wu_a = it.wu + alpha*del_wu_a;
            tl_a = it.tl + alpha*del_tl_a;
            tu_a = it.tu + alpha*del_tu_a;
            zl_a = it.zl + alpha*del_zl_a;
            zu_a = it.zu + alpha*del_zu_a;
            mu_aff = (sl_a'*wl_a + su_a'*wu_a + tl_a'*zl_a + tu_a'*zu_a)/(2*(n+r));
    
            % set the centering parameter sigma
            sigma = (mu_aff/mu)^3;
    
            % calculate new step
            psil = it.bound_cl*(it.zl - sigma*mu*help.Tl1*er + help.Tl1*diag(del_tl_a)*diag(del_zl_a)*er);
            psiu = it.bound_cu*(it.zu - sigma*mu*help.Tu1*er + help.Tu1*diag(del_tu_a)*diag(del_zu_a)*er);
            phil = it.bound_xl*(it.wl - help.Sl1*sigma*mu*en + help.Sl1*diag(del_sl_a)*diag(del_wl_a)*en);
            phiu = it.bound_xu*(it.wu - help.Su1*sigma*mu*en + help.Su1*diag(del_su_a)*diag(del_wu_a)*en);
            
            omega_1 = -nlp.H*it.x -nlp.c + nlp.A'*it.y + nlp.C'*(it.zl-it.zu) + it.bound_xl*it.wl - it.bound_xu*it.wu ...
                - it.bound_xl*(phil + help.Sl1*help.Wl*help.beta_l) + it.bound_xu*(phiu + help.Su1*help.Wu*help.beta_u);
            xi = -psil + psiu - it.bound_cl*help.Tl1*help.Zl*help.rho_l + it.bound_cu*help.Tu1*help.Zu*help.rho_u;
            omega_3 = help.PSI1*xi;
            omega = [omega_1;omega_2;omega_3];
            sol = mat\omega;
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
            
            % calculate step length alpha
            step = [del_sl; del_su; del_wl; del_wu; del_tl; del_tu; del_zl; del_zu];
            curr = [it.sl; it.su; it.wl; it.wu; it.tl; it.tu; it.zl; it.zu];
            [alpha] = step_length(step, curr, eta);
            
            % update iterate
            it.x = it.x + alpha*del_x;
            it.y = it.y + alpha*del_y;
            it.zl = it.zl + alpha*del_zl;
            it.zu = it.zu + alpha*del_zu;
            it.wl = it.wl + alpha*del_wl;
            it.wu = it.wu + alpha*del_wu;
            it.sl = it.sl + alpha*del_sl;
            it.su = it.su + alpha*del_su;
            it.tl = it.tl + alpha*del_tl;
            it.tu = it.tu + alpha*del_tu;

            mu_n = (it.sl'*it.wl + it.su'*it.wu + it.tl'*it.zl + it.tu'*it.zu)/(2*(n+r));
            obj = cutest_obj(it.x);

            data = [data; iterations obj mu_n mu-mu_n];
        if mu_n < 1*10^(-15)
            break;
        end
    end
    fprintf('      Iteration  objective     mu-value  difference in mu\n');
    disp(data);
    cutest_terminate
end