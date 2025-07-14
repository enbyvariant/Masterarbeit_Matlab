% Implementation of the Interior Points Method for LPs

function [x,sl,su,tl,tu,y,wl,wu,zl,zu, obj, iterations] = Interior_Points_LP(iter, A, C, cl, cu, xl, xu, b, c)
    
    % Initial values
    n = size(A, 2);
    m = size(A, 1);
    r = size(C, 1);
    eta = .9;
    en = ones(n,1);
    er = ones(r,1);

    % Compute starting point
    M = [A zeros(m,2*n+2*r);
        eye(n) -eye(n) zeros(n,n+2*r);
        -eye(n) zeros(n) -eye(n) zeros(n,2*r);
        C zeros(r,2*n) -eye(r) zeros(r);
        -C zeros(2*n+r) -eye(r)];
    vector = [b;cl;-cu;xl;-xu];
    a_pri = M' * linsolve(M*M',vector);
    x = a_pri(1:n);
    sl = a_pri(n+1:2*n);
    su = a_pri(2*n+1:3*n);
    tl = a_pri(3*n+1:3*n+r);
    tu = a_pri(3*n+r+1:3*n+2*r);

    N = [A C -C];
    yz = linsolve(N*N',N*c);
    y = yz(1:m);
    zl = yz(m+1:m+r);
    zu = yz(m+r+1:m+2*r);
    wlu = c - N'*yz;
    wl = 1/2*wlu;
    wu = -1/2*wlu;
    
    delta_pri = max([-3/2 *sl; -3/2*su; -3/2*tl;-3/2*tu;0]);
    delta_dual = max([-3/2 *wl; -3/2*wu; -3/2*zl; -3/2*zu; 0]);

    sl = sl + delta_pri * en;
    su = su + delta_pri * en;
    tl = tl + delta_pri * er;
    tu = tu + delta_pri * er;
    wl = wl + delta_dual * en;
    wu = wu + delta_dual * en;
    zl = zl + delta_dual * er;
    zu = zu + delta_dual * er;
    
    delta_pri = 1/2*(sl' * wl + su'*wu + tl'*zl + tu'*zu)/(en' * wl + en'*wu + er'*zl + er'*zu);
    delta_dual = 1/2*(sl' * wl + su'*wu + tl'*zl + tu'*zu)/(en' * sl + en'*su + er'*tl + er'*tu);

    sl = sl + delta_pri * en;
    su = su + delta_pri * en;
    tl = tl + delta_pri * er;
    tu = tu + delta_pri * er;
    wl = wl + delta_dual * en;
    wu = wu + delta_dual * en;
    zl = zl + delta_dual * er;
    zu = zu + delta_dual * er;

    % Main Loop
    for iterations = 1:iter

        % Update variables
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
        Wl = diag(wl);
        Wu = diag(wu);
        Zl = diag(zl);
        Zu = diag(zu);

        PHI = Sl1*Wl + Su1*Wu;
        
        r_c = A' * y + s - c;
        r_b = A * x - b;
        mat = [zeros(n,n) A' eye(n);A zeros(m,m) zeros(m,n); S zeros(n,m) X];

        % Solve for the affine search direction
         sol = linsolve(mat, [-r_c;-r_b;-X*S*en]);

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
        sol = linsolve(mat, [-r_c;-r_b; -X * S * en - diag(delta_x) * diag(delta_s) * en + sigma * mu * en]);
        
         delta_x = sol(1:n);
         delta_lambda = sol(n+1:m+n);
         delta_s = sol(m+n+1:2*n+m);

        % Determine primal and dual step length
        alpha_pri = min([1 eta * alpha_pri_max]);
        alpha_dual = min([1 eta * alpha_dual_max]);

        eta = eta + (1-eta)/2;

        % Calculate next iterate
        x = x + alpha_pri * delta_x;
        y = y + alpha_dual * delta_lambda;
        s = s + alpha_dual * delta_s;
        
        % Convergence
        if mu < .00001
            break
        end
    end
end