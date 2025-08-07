function [x,y,wl, wu, sl,su,tl,tu,zl,zu, mu_n, obj, iterations, data] = Interior_Points_LP_alt(iter, A, b, xl, xu, c, cl, cu, C, dim)
    
    % Initial values
    n = dim.n;
    m = dim.m;
    r = dim.r;
    en = ones(n,1);
    er = ones(r,1);
    eta = 0.995;

    % Compute starting point
    [x, sl, su, tl, tu, y, wl, wu, zl, zu, bound_xl, bound_xu, bound_cl, bound_cu] = Init_LP(A, C, cl, cu, xl, xu, b, c, dim);
    iterations = 0;
    data = zeros(iter,4);
    p = cutest_setup;


    while 1
        if iterations > iter
            break;
        end
        iterations = iterations + 1;

            % update helping variables
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

            rhol = bound_cl*(C*x - cl -tl);
            rhou = bound_cu*(-C*x + cu -tu);
            betal = bound_xl*(x -xl -sl);
            betau = bound_xu*(-x + xu -su);
            
            %set up the affine equation system
            phi = bound_xl*Sl1*Wl + bound_xu*Su1*Wu;
            psi = bound_cl*Tl1*Zl + bound_cu*Tu1*Zu;
            psi1 = zeros(r);
            for i = 1:r
                if psi(i,i) ~= 0
                psi1(i,i) = 1/psi(i,i);
                end
            end
            mat = [phi     A'         C';
                    A    zeros(m)   zeros(m,r);
                    C    zeros(r,m)   -psi1   ];
            omega_1 = - c + A'*y + C'*(zl - zu) - bound_xl*Sl1*Wl*betal + bound_xu*Su1*Wu*betau;
            omega_2 = -(A*x -b);
            xi = -zl - bound_cl*Tl1*Zl*rhol + zu + bound_cu*Tu1*Zu*rhou;
            omega_3 = psi1*xi;
            omega = [omega_1; omega_2; omega_3];

            % calculate affine step
            sol = mat\omega;
            del_x_a = sol(1:n);
    
            del_sl_a = bound_xl*(del_x_a + betal);
            del_su_a = bound_xu*(-del_x_a + betau);
            del_wl_a = bound_xl*(-wl -Sl1*Wl*del_sl_a);
            del_wu_a = bound_xu*(-wu - Su1*Wu*del_su_a);
        
            del_tl_a = bound_cl*(C*del_x_a + rhol);
            del_tu_a = bound_cu*(-C*del_x_a + rhou);
            del_zu_a = bound_cu*(-zu -Tu1*Zu*del_tu_a);
            del_zl_a = bound_cl*(-zl - Tl1*Zl*del_tl_a);
            
            % calculate duality measure
            mu = (sl'*wl + su'*wu + tl'*zl + tu'*zu)/(2*(n+r));
            
            % calculate primal affine step length

            step = [del_sl_a; del_su_a; del_wl_a; del_wu_a];
            curr = [sl; su; wl; wu];
            alpha_pri = step_length(step,curr,1);

            % calculate dual affine step length
            step = [del_tl_a; del_tu_a; del_zl_a; del_zu_a];
            curr = [tl; tu; zl; zu];
            alpha_dual = step_length(step, curr, eta);
            %------------------------
            % calculate affine duality measure
            sl_a = sl + alpha_pri*del_sl_a;
            su_a = su + alpha_pri*del_su_a;
            wl_a = wl + alpha_pri*del_wl_a;
            wu_a = wu + alpha_pri*del_wu_a;
            tl_a = tl + alpha_dual*del_tl_a;
            tu_a = tu + alpha_dual*del_tu_a;
            zl_a = zl + alpha_dual*del_zl_a;
            zu_a = zu + alpha_dual*del_zu_a;
            mu_aff = (sl_a'*wl_a + su_a'*wu_a + tl_a'*zl_a + tu_a'*zu_a)/(2*(n+r));
    
            % set the centering parameter sigma
            sigma = (mu_aff/mu)^3;
    
            % set up equation system
            psil = bound_cl*(zl - sigma*mu*Tl1*er + Tl1*diag(del_tl_a)*diag(del_zl_a)*er);
            psiu = bound_cu*(zu - sigma*mu*Tu1*er + Tu1*diag(del_tu_a)*diag(del_zu_a)*er);
            phil = bound_xl*(wl - Sl1*sigma*mu*en + Sl1*diag(del_sl_a)*diag(del_wl_a)*en);
            phiu = bound_xu*(wu - Su1*sigma*mu*en + Su1*diag(del_su_a)*diag(del_wu_a)*en);
            
            omega_1 = -c + A'*y + C'*(zl-zu) + bound_xl*wl - bound_xu*wu - bound_xl*(phil + Sl1*Wl*betal) + bound_xu*(phiu + Su1*Wu*betau);
            xi = -psil + psiu - bound_cl*Tl1*Zl*rhol + bound_cu*Tu1*Zu*rhou;
            omega_3 = psi1*xi;
            omega = [omega_1;omega_2;omega_3];

            % calculate new step
            sol = mat\omega;
            del_x = sol(1:n);
            del_y = -sol(n+1:n+m);
            del_sl = bound_xl*(del_x + betal);
            del_su = bound_xu*(-del_x + betau);
            del_wl = bound_xl*(-phil -Sl1*Wl*del_sl);
            del_wu = bound_xu*(-phiu - Su1*Wu*del_su);
        
            del_tl = bound_cl*(C*del_x + rhol);
            del_tu = bound_cu*(-C*del_x + rhou);
            del_zu = bound_cu*(-psiu -Tu1*Zu*del_tu);
            del_zl = bound_cl*(-psil - Tl1*Zl*del_tl);
            
            %----------------------------
            % calculate primal step length

            step = [del_sl; del_su; del_wl; del_wu];
            curr = [sl; su; wl; wu];
            alpha_pri = step_length(step,curr,eta);

            % calculate dual step length
            step = [del_tl; del_tu; del_zl; del_zu];
            curr = [tl; tu; zl; zu];
            alpha_dual = step_length(step, curr, eta);

            % update iterate
            x = x + alpha_pri*del_x;
            y = y + alpha_dual*del_y;
            zl = zl + alpha_dual*del_zl;
            zu = zu + alpha_dual*del_zu;
            wl = wl + alpha_pri*del_wl;
            wu = wu + alpha_pri*del_wu;
            sl = sl + alpha_pri*del_sl;
            su = su + alpha_pri*del_su;
            tl = tl + alpha_dual*del_tl;
            tu = tu + alpha_dual*del_tu;

            mu_n = (sl'*wl + su'*wu + tl'*zl + tu'*zu)/(2*(n+r));
            obj = cutest_obj(x);

            data(iterations, :) = [iterations obj mu_n mu-mu_n];
        if mu_n < 1*10^(-15)
            break;
        end
    end
    fprintf('      Iteration  objective     mu-value  difference in mu\n');
    disp(data);
    cutest_terminate
end