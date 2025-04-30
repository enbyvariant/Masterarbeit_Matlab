% Implementation of the Interior Points Method for LPs
% jetzt auch in github

function [x,y,v_1,v_2, mu] = Interior_Points_QP(iter, A, b, xl, xu, H, c, cl, cu, C)
    
    % Initial values
    n = size(A, 2);
    m = size(A, 1);
    r = size(C,1);
    e = ones(n,1);
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

        Hphi = H + Sl1*Wl + Su1*Wu;
        psi = Tl1*Zl + Tu1*Zu;
        rhol = (C*x - cl -tl)*bound_cl;
        rhou = (-C*x + cu -tu)*bound_cu;
        betal = (x -xl -sl)*bound_xl;
        betau = (-x + xu -su)*bound_xu;
        xi = zl - zu - Tu1*Zu*(rhol+rhou);
        foo = xi\psi;
        mat = [Hphi A' C';
            A zeros(m,m+r);
            C zeros(r,m) inv(psi)];
        omega = [H*x + c - A'*y - bound_cl*C'*zl + bound_cu*C'*zu +Sl1*Wl*betal -Su1*Wu*betau;
            A*x-b;
            rhol + foo];
        sol = omega\mat;
        del_x = sol(1:n);
        del_z = sol(n+m+1:n+m+r);
        del_sl = bound_xl*(del_x + betal);
        del_su = bound_xu*(-del_x + betau);
        del_wl = bound_xl*(-wl -Sl1*Wl*del_sl);
        del_wu = bound_xu*(-wu - Su1*Wu*del_su);
    
        del_tl = (-del_z -xi)\psi;
        del_tu = -del_tl + rhou + rhol;
        del_zu = -zu -Tu1*Zu*del_tu;
        del_zl = del_z + del_zu;
        
        % calculate duality measure
        mu = (sl'*wl + su'*wu + tl'*zl + tu'*zu)/(2*n+2*r);
        
        % calculate affine step length
        step = [del_sl; del_su; del_wl; del_wu; del_tl; del_tu; del_zl; del_zu];
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
        mu_aff = ((sl+ alpha*del_sl)' * (wl + alpha*del_wl) + (su + alpha*del_su)' * (wu + del_wu) + (tl + alpha*del_tl)'*(zl + alpha*del_zl) + (tu + alpha*del_tu)'*(zu + alpha*del_zu))/(2*n + 2*r);

        % set the centering parameter sigma
        sigma = (mu_aff/mu)^3;
    end



    switch bound
        case "lu"
            for i = 1:iter

                Sl1 = diag(ones(n,1)./sl);
                Wl = diag(wl);
                Su1 = diag(ones(n,1)./su);
                Wu = diag(wu);
                Hphi = H + Sl1*Wl + Su1*Wu;
    
                %compute duality measure
                mu = (sl' * wl + su' * wu)/(2*n);
    
                %compute affine step direction
                mat = [Hphi A'; A zeros(m)];
                omega = [H*x + c -A'*y - wl + wu + Sl1*Wl*(x -xl -sl) -Su1*Wu*(-x +xu -su) + wl - wu; A*x - b];
                sol = linsolve(mat, -omega);
                delta_x = sol(1:n);
                delta_sl = delta_x + (x -xl -sl);
                delta_su = -delta_x + (-x +xu -su);
                delta_wl = -Sl1*Wl*delta_x - Sl1*Wl*(x -xl -sl) - wl;
                delta_wu = Su1*Wu*delta_x - Su1*Wu*(-x +xu -su) - wu;
                
                %compute affine step length alpha_aff
                step = [delta_sl; delta_su; delta_wl; delta_wu];
                curr = [sl; su; wl; wu];
                index = find(step < 0);
                if isempty(index)
                    alpha = 1;
                else
                    alpha = min(curr(index)./(-step(index)));
                    if alpha > 1
                        alpha = 1;
                    end
                end
    
                % compute affine duality measure
                mu_aff = ((sl'+ alpha*delta_sl) * (wl + alpha*delta_wl) + (su' + alpha*delta_su) * (wu + delta_wu))/(2*n);
                
                % set the centering parameter sigma
                sigma = (mu_aff/mu)^3;

                % compute new step direction
                omega = [H*x + c -A'*y - wl + wu + Sl1*Wl*(x -xl -sl) -Su1*Wu*(-x +xu -su) + (wl - Sl1*sigma*mu*e + Sl1*diag(delta_sl)*diag(delta_wl)) - (wu - Sl1*sigma*mu*e + Su1*diag(delta_su)*diag(delta_wu)); A*x - b];
                sol = linsolve(mat,-omega);
                delta_x = sol(1:n);
                delta_y = sol(n+1:n+m);
                delta_sl = delta_x + (x -xl -sl);
                delta_su = -delta_x + (-x +xu -su);
                delta_wl = -Sl1*Wl*delta_x - Sl1*Wl*(x -xl -sl) - wl;
                delta_wu = Su1*Wu*delta_x - Su1*Wu*(-x +xu -su) - wu;

                % compute new step length alpha
                step = [delta_sl; delta_su; delta_wl; delta_wu];
                curr = [sl; su; wl; wu];
                index = find(step < 0);
                if isempty(index)
                    alpha = 1;
                else
                    alpha = min(-tao*curr(index)./(step(index)));
                    if alpha > 1
                        alpha = 1;
                    end
                end
                
                % update iterate
                x = x + alpha*delta_x;
                y = y + alpha*delta_y;
                sl = sl + alpha*delta_sl;
                su = su + alpha*delta_su;
                wl = wl + alpha*delta_wl;
                wu = wu + alpha*delta_wu;
            end
        case "l"

        case "u"

        otherwise
            for i = 1:iter
                
            end
    end


end