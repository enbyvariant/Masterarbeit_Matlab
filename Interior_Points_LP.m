% Implementation of the Interior Points Method for LPs

function [x,lambda,s] = Interior_Points_LP(iter, A, b, c)
    
    % Initial values
    n = size(A, 2);
    m = size(A, 1);
    eta = .9;
    e = ones(n,1);

    % Compute starting point
    x = A' * linsolve(A*A',b);
    lambda = linsolve(A*A', A*c);
    s = c - A'*lambda;
    
    del_x = max([-3/2 * min(x) 0]);
    del_s = max([-3/2 * min(s) 0]);

    x = x + del_x * e;
    s = s + del_s * e;
    
    del_x = (x' * s)/(e' * s)/2;
    del_s = (x' * s)/(e' * x)/2;

    x = x + del_x * e;
    s = s + del_s * e;

    % Main Loop
    for i = 1:iter

        % Update variables
        X = diag(x);
        S = diag(s);
        r_c = A' * lambda + s - c;
        r_b = A * x - b;
        mat = [zeros(n,n) A' eye(n);A zeros(m,m) zeros(m,n); S zeros(n,m) X];

        % Solve for the affine search direction
         sol = linsolve(mat, [-r_c;-r_b;-X*S*e]);

         delta_x = sol(1:n);
         %delta_lambda = sol(n+1:m+n);
         delta_s = sol(m+n+1:2*n+m);
        
        % Determine affine primal step length alpha_pri_aff
        alpha_pri_max = Inf;
        for j = 1:n
            if delta_x(j) < 0
                if -x(j)/delta_x(j) < alpha_pri_max
                    alpha_pri_max = -x(j)/delta_x(j);
                end
            end
        end
        alpha_pri = min([1 alpha_pri_max]);
        
        % Determine affine dual step length alpha_dual_aff
        alpha_dual_max = Inf;
        for j = 1:n
            if delta_s(j) < 0
                if -s(j)/delta_s(j) < alpha_dual_max
                    alpha_dual_max = -s(j)/delta_s(j);
                end
            end
        end
        alpha_dual = min([1 alpha_dual_max]);
        
        % Calculate duality measure mu
        mu = x' * s /n;

        % Calculate affine duality measure mu_aff
        mu_aff = (x + alpha_pri * delta_x)' * (s + alpha_dual * delta_s)/n;
        
        % Set the centering parameter
        sigma = (mu_aff/mu)^3;
        
        % Calculate the search direction 
        sol = linsolve(mat, [-r_c;-r_b; -X * S * e - diag(delta_x) * diag(delta_s) * e + sigma * mu * e]);
        
         delta_x = sol(1:n);
         delta_lambda = sol(n+1:m+n);
         delta_s = sol(m+n+1:2*n+m);

        % Determine primal and dual step length
        alpha_pri = min([1 eta * alpha_pri_max]);
        alpha_dual = min([1 eta * alpha_dual_max]);

        eta = eta + (1-eta)/2;

        % Calculate next iterate
        x = x + alpha_pri * delta_x;
        lambda = lambda + alpha_dual * delta_lambda;
        s = s + alpha_dual * delta_s;
        
        % Convergence
        if mu < .00001
            break
        end
    end
end