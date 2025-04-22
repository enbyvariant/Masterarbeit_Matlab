% Implementation of the Interior Points Method for LPs
% jetzt auch in github

function [x,y,v_1,v_2, mu] = Interior_Points_QP(iter, A, b, xl, xu, H, c, x, y)
    
    % Initial values
    n = size(A, 2);
    m = size(A, 1);
    e = ones(n,1);
    tao = .5;

    % Compute starting point
    [x, sl, su, y, wl, wu, bound] = Interior_Points_Init(H, A, b, c, xl, xu, x, y);
    

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

    end


end