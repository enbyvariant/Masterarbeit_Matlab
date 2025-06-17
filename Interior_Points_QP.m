% Implementation of the Interior Points Method for LPs
% jetzt auch in github

function [x,y,wl, wu, sl,su, mu_n, obj, iterations] = Interior_Points_QP(iter, A, b, xl, xu, H, c, cl, cu, C)
    
    % Initial values
    n = size(A, 2);
    m = size(A, 1);
    r = size(C,1);
    en = ones(n,1);
    er = ones(r,1);
    tao = 0.9;

    % Compute starting point
    [x, sl, su, tl, tu, y, wl, wu, zl, zu, bound_xl,bound_xu, bound_cl, bound_cu] = Interior_Points_Init(H, A, b, c, xl, xu, cl, cu, C);
    % x = [1;1;1;1];
    % sl = x-xl;
    % su = zeros(4,1);
    % y = ones(2,1);
    % wl = ones(4,1);
    % wu = ones(4,1);
    % bound_xl = eye(4);
    % bound_xu = zeros(4);
    iterations = 0;

    data = [];
    p = cutest_setup;


    while 1
        if iterations > iter
            break;
        end
        iterations = iterations + 1;

        % if inequalities do not exist, use smaller equation system
        if r == 0
           % Case: Inequalities do not exist
           %----------------------------------

           % update  helping variables
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
            betal = bound_xl*(x -xl -sl);
            betau = bound_xu*(-x + xu -su);

            Hphi = H + bound_xl*Sl1*Wl + bound_xu*Su1*Wu;
            mat = [Hphi     A';
                    A    zeros(m)];
            omega_1 = -H*x - c + A'*y - bound_xl*Sl1*Wl*betal + bound_xu*Su1*Wu*betau;
            omega_2 = -(A*x -b);
            omega = [omega_1; omega_2];
            sol = mat\omega;
            del_x_a = sol(1:n);
            %del_y_a = -sol(n+1:n+m);

            % calculate other variables of expanded system
            del_sl_a = bound_xl*(del_x_a + betal);
            del_su_a = bound_xu*(-del_x_a + betau);
            del_wl_a = bound_xl*(-wl -Sl1*Wl*del_sl_a);
            del_wu_a = bound_xu*(-wu - Su1*Wu*del_su_a);

            % calculate affine step length
            step = [del_sl_a; del_su_a; del_wl_a; del_wu_a];
            curr = [sl; su; wl; wu];
            index = find(step < 0);
            if isempty(index)
                alpha_aff = 1;
            else
                alpha_aff = min(-tao*curr(index)./step(index));
                if alpha_aff > 1
                    alpha_aff = 1;
                end
            end
            
            % calculate provisional new step
            %x_a = x + alpha_aff*del_x_a;
            %y_a = y + alpha_aff*del_y_a;
            sl_a = bound_xl*(sl + alpha_aff*del_sl_a);
            su_a = bound_xu*(su + alpha_aff*del_su_a);
            wl_a = bound_xl*(wl + alpha_aff*del_wl_a);
            wu_a = bound_xu*(wu + alpha_aff*del_wu_a);

            % calculate duality measure mu
            mu = (sl'*wl + su'*wu)/(2*n);

            % calculate affine duality measure mu_aff
            mu_aff = (sl_a'*wl_a + su_a'*wu_a)/(2*n);

            % calculate centering parameter sigma
            sigma = (mu_aff/mu)^3;
            if sigma > 1
                sigma = 1;
            end
            
            phil = wl - Sl1*sigma*mu*en + Sl1*diag(del_sl_a)*diag(del_wl_a)*en;
            phiu = wu - Su1*sigma*mu*en + Su1*diag(del_su_a)*diag(del_wu_a)*en;
            omega_1 = -H*x -c + A'*y + bound_xl*wl - bound_xu*wu - bound_xl*(phil + Sl1*Wl*betal) + bound_xu*(phiu + Su1*Wu*betau);
            omega = [omega_1; -(A*x - b)];
            sol = mat\omega;
            del_x = sol(1:n);
            del_y = -sol(n+1:n+m);

            % calculate all variables of the expanded system
            del_sl = bound_xl*(del_x + betal);
            del_su = bound_xu*(-del_x + betau);
            del_wl = bound_xl*(-Sl1*Wl*del_sl - phil);
            del_wu = bound_xu*(-Su1*Wu*del_su - phiu);

            % calculate step length alpha
            step = [del_sl; del_su; del_wl; del_wu];
            curr = [sl; su; wl; wu];
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

            sl = bound_xl*(sl + alpha*del_sl);
            su = bound_xu*(su + alpha*del_su);

            wl = bound_xl*(wl + alpha*del_wl);
            wu = bound_xu*(wu + alpha*del_wu);
            

            mu_n = (sl'*wl+su'*wu)/(2*n);
            obj = cutest_obj(x);
            data = [data; iterations obj mu_n mu-mu_n];

        else
            % Case Inequalities do exist
            %----------------------------%

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
            
            % calculate affine step
            Hphi = H + bound_xl*Sl1*Wl + bound_xu*Su1*Wu;
            psi = bound_cl*Tl1*Zl + bound_cu*Tu1*Zu;
            psi1 = zeros(r);
            for i = 1:r
                psi1(i,i) = 1/psi(i,i);
            end
            mat = [Hphi     A'         C';
                    A    zeros(m)   zeros(m,r);
                    C    zeros(r,m)   -psi1   ];
            omega_1 = -H*x - c + A'*y - bound_xl*Sl1*Wl*betal + bound_xu*Su1*Wu*betau;
            omega_2 = -(A*x -b);
            xi = -zl - bound_cl*Tl1*Zl*rhol + zu + bound_cu*Tu1*Zu*rhou;
            omega_3 = psi1*xi;
            omega = [omega_1; omega_2; omega_3];
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
            sl_a = sl + alpha*del_sl_a;
            su_a = su + alpha*del_su_a;
            wl_a = wl + alpha*del_wl_a;
            wu_a = wu + alpha*del_wu_a;
            tl_a = tl + alpha*del_tl_a;
            tu_a = tu + alpha*del_tu_a;
            zl_a = zl + alpha*del_zl_a;
            zu_a = zu + alpha*del_zu_a;
            mu_aff = (sl_a'*wl_a + su_a'*wu_a + tl_a'*zl_a + tu_a'*zu_a)/(2*(n+r));
    
            % set the centering parameter sigma
            sigma = (mu_aff/mu)^3;
    
            % calculate new step
            psil = bound_cl*(zl - sigma*mu*Tl1*er + Tl1*diag(del_tl_a)*diag(del_zl_a)*er);
            psiu = bound_cu*(zu - sigma*mu*Tu1*er + Tu1*diag(del_tu_a)*diag(del_zu_a)*er);
            phil = bound_xl*(wl - Sl1*sigma*mu*en + Sl1*diag(del_sl_a)*diag(del_wl_a)*en);
            phiu = bound_xu*(wu - Su1*sigma*mu*en + Su1*diag(del_su_a)*diag(del_wu_a)*en);
            omega_1 = -H*x -c + A'*y + bound_xl*wl - bound_xu*wu - bound_xl*(phil + Sl1*Wl*betal) + bound_xu*(phiu + Su1*Wu*betau);
            xi = -psil + psiu - bound_cl*Tl1*Zl*rhol + bound_cu*Tu1*Zu*rhou;
            omega_3 = psi1*xi;
            omega = [omega_1;omega_2;omega_3];
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

            mu_n = (sl'*wl + su'*wu + tl'*zl + tu'*zu)/(2*(n+r));
            obj = cutest_obj(x);

            data = [data; iterations obj mu_n mu-mu_n];
        end
        if mu_n < 1*10^(-14)
            break;
        end
    end
    fprintf('      Iteration  objective    mu-value  difference in mu\n');
    disp(data);
    cutest_terminate
end