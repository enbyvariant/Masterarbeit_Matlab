% Implementation of the Interior Points Method for LPs
% jetzt auch in github

function [x,y,v_1,v_2, mu] = Interior_Points_QP(iter, A, b, xl, xu, H, c, x, y, v_1,v_2)
    
    % Initial values
    n = size(A, 2);
    m = size(A, 1);
    e = ones(n,1);
    tao = .5;

    % Compute starting point
    V_1 = diag(v_1);
    V_2 = diag(v_2);
    Xl = diag(xl);
    Xu = diag(xu);
    X = diag(x);
    mat = [H -A' -eye(n) eye(n); A zeros(m,m) zeros(m,2*n); V_1 zeros(n,m) X-Xl zeros(n,n); -V_2 zeros(n,m+n) -X+Xu];
    sol = linsolve(mat, [-H*x + A'*y + v_1-v_2 - c; -A*x + b; -V_1*(X-Xl)*e; V_2*(X-Xu)*e]);
    delta_x = sol(1:n);
    %delta_y = sol(n+1:m+n);
    delta_v_1 = sol(m+n+1:2*n+m);
    delta_v_2 = sol(m+2*n+1:3*n+m);
    x = max(ones(n,1),abs(x+delta_x));
    v_1 = max(ones(n,1),abs(v_1+delta_v_1));
    v_2 = max(ones(n,1),abs(v_2+delta_v_2));

    % Main Loop
    for i = 1:iter

        % Update variables
        X = diag(x);
        V_1 = diag(v_1);
        V_2 = diag(v_2);
        r_d = H*x - A' * y - v_1 + v_2 + c;
        r_p = A * x - b;
        mat = [H -A' -eye(n) eye(n); A zeros(m,m) zeros(m,2*n); V_1 zeros(n,m) X-Xl zeros(n); -V_2 zeros(n,m+n) -X+Xu];

        % Solve for the affine search direction
        sol = linsolve(mat, [-r_d;-r_p;-V_1*(X-Xl)*e; V_2*(X-Xu)*e]);

        delta_x = sol(1:n);
        %delta_lambda = sol(n+1:m+n);
        delta_v_1 = sol(m+n+1:2*n+m);
        delta_v_2 = sol(m+2*n+1:3*n+m);
        
        % Determine affine step length alpha_aff
        del = [delta_v_1; delta_v_2];
        z = [v_1;v_2];
        index = find(del < 0);
        if isempty(index)
            alpha = 1;
        else
            alpha = min(z(index)./(-del(index)));
            if alpha > 1
            alpha = 1;
            end
        end
        
        % Calculate duality measure mu
        mu_1 = (x-xl)' * v_1 /n;
        mu_2 = (-x+xu)' * v_2 /n;
        mu = (mu_1 + mu_2)/2;

        % Calculate affine duality measure mu_aff
        mu_aff_1 = (x + alpha*delta_x-xl)' * (v_1 + alpha * delta_v_1)/n;
        mu_aff_2 = (-x - alpha*delta_x+xu)' * (v_2 + alpha * delta_v_2)/n;
        
        % Set the centering parameter
        sigma_1 = (mu_aff_1/mu_1)^3;
        sigma_2 = (mu_aff_2/mu_2)^3;

        % Calculate the search direction 
        sol = linsolve(mat, [-r_d;-r_p; -V_1*(X-Xl)*e - diag(delta_x) * diag(delta_v_1) * e + sigma_1 * mu_1 * e; V_2*(X-Xu)*e + diag(delta_x) * diag(delta_v_2) * e + sigma_2 * mu_2 * e]);
        
        delta_x = sol(1:n);
        delta_y = sol(n+1:m+n);
        delta_v = sol(m+n+1:3*n+m);

        % Determine step length
        % index = find(delta_x < 0);
        % if isempty(index)
        %     alpha_pri = 1;
        % else
        %     alpha_pri = min(-tao*x(index)./delta_x(index));
        %     if alpha_pri > 1
        %         alpha_pri = 1;
        %     end
        % end
        v = [v_1; v_2];
        index = find(delta_v < 0);
        if isempty(index)
            alpha_dual = 1;
        else
            alpha_dual = min(-tao*v(index)./delta_v(index));
            if alpha_dual > 1
                alpha_dual = 1;
            end
        end
        alpha = alpha_dual; %min([alpha_pri alpha_dual]);

        tao = tao + (1-tao)/2;

        % Calculate next iterate
        x = x + alpha * delta_x;
        y = y + alpha * delta_y;
        v_1 = v_1 + alpha * delta_v(1:n);
        v_2 = v_2 + alpha * delta_v(n+1:2*n);
        
        % Convergence
        if mu < .00000000000000000000001
            break
        end
    end

end