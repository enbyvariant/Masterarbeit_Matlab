function [it, mu_n, obj, iterations, data] = Interior_Points_LP_alt(iter, A, b, xl, xu, c, cl, cu, C, nlp,dim)
    
    % Initial values
    n = size(A, 2);
    m = size(A, 1);
    r = size(C,1);
    en = ones(n,1);
    er = ones(r,1);
    eta = 0.995;

    % Compute starting point
    [it] = Init_LP(nlp,dim);
    iterations = 0;
    data = [];
    p = cutest_setup;


    while 1
        if iterations > iter
            break;
        end
        iterations = iterations + 1;

            % update helping variables
            help = helpers(dim,nlp,it);
            % Sl1 = zeros(n);
            % Su1 = zeros(n);
            % for i = 1:n
            %     if sl(i) ~= 0
            %         Sl1(i,i) = 1/sl(i);
            %     end
            %     if su(i) ~= 0
            %         Su1(i,i) = 1/su(i);
            %     end
            % end
            % Wl = diag(wl);
            % Wu = diag(wu);
            % 
            % Tl1 = zeros(r);
            % Tu1 = zeros(r);
            % for i = 1:r
            %     if tl(i) ~= 0
            %         Tl1(i,i) = 1/tl(i);
            %     end
            %     if tu(i) ~= 0
            %         Tu1(i,i) = 1/tu(i);
            %     end
            % end
            % Zl = diag(zl);
            % Zu = diag(zu);
            % rhol = bound_cl*(C*x - cl -tl);
            % rhou = bound_cu*(-C*x + cu -tu);
            % betal = bound_xl*(x -xl -sl);
            % betau = bound_xu*(-x + xu -su);
            
            %set up the affine equation system
            phi = it.bound_xl*help.Sl1*help.Wl + it.bound_xu*help.Su1*help.Wu;
            psi = it.bound_cl*help.Tl1*help.Zl + it.bound_cu*help.Tu1*help.Zu;
            psi1 = zeros(r);
            for i = 1:r
                if psi(i,i) ~= 0
                psi1(i,i) = 1/psi(i,i);
                end
            end
            mat = [phi     A'         C';
                    A    zeros(m)   zeros(m,r);
                    C    zeros(r,m)   -psi1   ];
            omega_1 = - c + A'*it.y + C'*(it.zl - it.zu) - it.bound_xl*help.Sl1*help.Wl*help.beta_l + it.bound_xu*help.Su1*help.Wu*help.beta_u;
            omega_2 = -(A*it.x -b);
            xi = -it.zl - it.bound_cl*help.Tl1*help.Zl*help.rho_l + it.zu + it.bound_cu*help.Tu1*help.Zu*help.rho_u;
            omega_3 = psi1*xi;
            omega = [omega_1; omega_2; omega_3];

            % calculate affine step
            sol = mat\omega;
            del_x_a = sol(1:n);
    
            del_sl_a = it.bound_xl*(del_x_a + help.beta_l);
            del_su_a = it.bound_xu*(-del_x_a + help.beta_u);
            del_wl_a = it.bound_xl*(-it.wl -help.Sl1*help.Wl*del_sl_a);
            del_wu_a = it.bound_xu*(-it.wu - help.Su1*help.Wu*del_su_a);
        
            del_tl_a = it.bound_cl*(C*del_x_a + help.rho_l);
            del_tu_a = it.bound_cu*(-C*del_x_a + help.rho_u);
            del_zu_a = it.bound_cu*(-it.zu -help.Tu1*help.Zu*del_tu_a);
            del_zl_a = it.bound_cl*(-it.zl - help.Tl1*help.Zl*del_tl_a);
            
            % calculate duality measure
            mu = (it.sl'*it.wl + it.su'*it.wu + it.tl'*it.zl + it.tu'*it.zu)/(2*(n+r));
            
            % calculate primal affine step length

            step = [del_sl_a; del_su_a; del_tl_a; del_tu_a];
            curr = [it.sl; it.su; it.tl; it.tu];

            index = find(step < 0);
            if isempty(index)
               alpha_pri = 1;
            else
                alpha_pri = min(curr(index)./(-step(index)));
                if alpha_pri > 1
                    alpha_pri = 1;
                end
            end

            % calculate dual affine step length
            step = [del_wl_a; del_wu_a; del_zl_a; del_zu_a];
            curr = [it.wl; it.wu; it.zl; it.zu];
            index = find(step < 0);
            if isempty(index)
               alpha_dual = 1;
            else
                alpha_dual = min(curr(index)./(-step(index)));
                if alpha_dual > 1
                    alpha_dual = 1;
                end
            end
            
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
            phil = it.bound_xl*(it.wl - help.Sl1*sigma*mu*en + help.Sl1*diag(del_sl_a)*diag(del_wl_a)*en);
            phiu = it.bound_xu*(it.wu - help.Su1*sigma*mu*en + help.Su1*diag(del_su_a)*diag(del_wu_a)*en);
            
            omega_1 = -c + A'*it.y + C'*(it.zl-it.zu) + it.bound_xl*it.wl - it.bound_xu*it.wu - it.bound_xl*(phil + help.Sl1*help.Wl*help.beta_l) + it.bound_xu*(phiu + help.Su1*help.Wu*help.beta_u);
            xi = -psil + psiu - it.bound_cl*help.Tl1*help.Zl*help.rho_l + it.bound_cu*help.Tu1*help.Zu*help.rho_u;
            omega_3 = psi1*xi;
            omega = [omega_1;omega_2;omega_3];

            % calculate new step
            sol = mat\omega;
            del_x = sol(1:n);
            del_y = -sol(n+1:n+m);
            del_sl = it.bound_xl*(del_x + help.beta_l);
            del_su = it.bound_xu*(-del_x + help.beta_u);
            del_wl = it.bound_xl*(-phil -help.Sl1*help.Wl*del_sl);
            del_wu = it.bound_xu*(-phiu - help.Su1*help.Wu*del_su);
        
            del_tl = it.bound_cl*(C*del_x + help.rho_l);
            del_tu = it.bound_cu*(-C*del_x + help.rho_u);
            del_zu = it.bound_cu*(-psiu -help.Tu1*help.Zu*del_tu);
            del_zl = it.bound_cl*(-psil - help.Tl1*help.Zl*del_tl);
            
            %----------------------------
            % calculate primal affine step length

            step = [del_sl; del_su; del_tl; del_tu];
            curr = [it.sl; it.su; it.tl; it.tu];

            index = find(step < 0);
            if isempty(index)
               alpha_pri = 1;
            else
                alpha_pri = eta*min(curr(index)./(-step(index)));
                if alpha_pri > 1
                    alpha_pri = 1;
                end
            end

            % calculate dual affine step length
            step = [del_wl; del_wu; del_zl; del_zu];
            curr = [it.wl; it.wu; it.zl; it.zu];
            index = find(step < 0);
            if isempty(index)
               alpha_dual = 1;
            else
                alpha_dual = eta*min(curr(index)./(-step(index)));
                if alpha_dual > 1
                    alpha_dual = 1;
                end
            end
            
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