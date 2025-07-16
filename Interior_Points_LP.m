
% Implementation of the Interior Points Method for LPs

function [x,sl,su,tl,tu,y,wl,wu,zl,zu,mu, iterations] = Interior_Points_LP(iter, A, C, cl, cu, xl, xu, b, c)
    
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
        -C zeros(r,2*n+r) -eye(r)];
    vector = [b;xl;-xu;cl;-cu];
    a_pri = M' * linsolve(M*M',vector);
    x = a_pri(1:n);
    sl = a_pri(n+1:2*n);
    su = a_pri(2*n+1:3*n);
    tl = a_pri(3*n+1:3*n+r);
    tu = a_pri(3*n+r+1:3*n+2*r);

    N = [A; C; -C];
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
    % x = [0;0];
    % sl = [1;1];
    % su = [1;1];
    % tl = [1;1];
    % tu = [1;1];
    % wl = [1;1];
    % wu = [1;1];
    % zl = [1;1];
    % zu = [1;1];
    % y = zeros(0,1);

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
        
        % Define helping variables
        PHI = Sl1*Wl + Su1*Wu;
        PSI = Tl1*Zl + Tu1*Zu;
        PSI1 = 1./PSI;

        beta_l = x - xl - sl;
        beta_u = -x + xu -su;
        rho_l = C*x - cl - tl;
        rho_u = -C*x + cu - tu;
        
        % set up equation system
        nu = -(c - A'*y - C'*(-zl + zu)) -Sl1*Wl*beta_l + Su1*Wu*beta_u;
        xi = PSI1*(-zl - Tl1*Zl*rho_l + zu + Tu1*Zu*rho_u);
        mat = [PHI A' C';
            A zeros(m,m) zeros(m,r);
            C zeros(r,m) -PSI1];

        % Solve for the affine search direction
         sol = linsolve(mat, [nu;-(A*x - b);xi]);

         delta_x = sol(1:n);
         delta_sl = delta_x + beta_l;
         delta_su = -delta_x + beta_u;
         delta_tl = C*delta_x + rho_l;
         delta_tu = -C*delta_x + rho_u;
         delta_wl = -Sl1*Wl*delta_sl -wl;
         delta_wu = -Su1*Wu*delta_su - wu;
         delta_zl = -Tl1*Zl*delta_tl - zl;
         delta_zu = -Tu1*Zu*delta_tu - zu;
        
        % Determine affine primal step length alpha_pri_aff
            step = [delta_sl; delta_su; delta_tl; delta_tu];
            curr = [sl; su; tl; tu];
            index = find(step < 0);
            if isempty(index)
                alpha_pri = 1;
            else
                alpha_pri = min(curr(index)./(-step(index)));
                if alpha_pri > 1
                    alpha_pri = 1;
                end
            end
        
        % Determine affine dual step length alpha_dual_aff
            step = [delta_zl; delta_zu; delta_wl; delta_wu];
            curr = [zl; zu; wl; wu];
            index = find(step < 0);
            if isempty(index)
                alpha_dual = 1;
            else
                alpha_dual = min(curr(index)./(-step(index)));
                if alpha_dual > 1
                    alpha_dual = 1;
                end
            end
        
        % Calculate duality measure mu
        mu = (sl'*wl + su'*wu + tl'*zl + tu'*zu)/(2*n+2*r);

        % Calculate affine duality measure mu_aff
        sl_a = sl + alpha_pri*delta_sl;
        su_a = su + alpha_pri*delta_su;
        tl_a = tl + alpha_pri*delta_tl;
        tu_a = tu + alpha_pri*delta_tu;
        wl_a = wl + alpha_dual*delta_wl;
        wu_a = wu + alpha_dual*delta_wu;
        zl_a = zl + alpha_dual*delta_zl;
        zu_a = zu + alpha_dual*delta_zu;
        mu_aff = (sl_a'*wl_a + su_a'*wu_a + tl_a'*zl_a + tu_a'*zu_a)/(2*n+2*r);
        
        % Set the centering parameter
        sigma = (mu_aff/mu)^3;
        
        % Calculate the search direction
        phi_l = wl - Sl1*(sigma*mu*en - diag(delta_sl)*delta_wl);
        phi_u = wu - Su1*(sigma*mu*en - diag(delta_su)*delta_wu);
        psi_l = zl - Tl1*(sigma*mu*er - diag(delta_tl)*delta_zl);
        psi_u = zu - Tu1*(sigma*mu*er - diag(delta_tu)*delta_zu);
        nu = -(c - A'*y - C'*(-zl + zu)- wl + wu) - phi_l - Sl1*Wl*beta_l + phi_u + Su1*Wu*beta_u;
        xi = PSI1*(-psi_l - Tl1*Zl*rho_l + psi_u + Tu1*Zu*rho_u);
        sol = linsolve(mat, [nu;-(A*x - b);xi]);
        
         delta_x = sol(1:n);
         delta_y = -sol(n+1:m+n);
         delta_sl = delta_x + beta_l;
         delta_su = -delta_x + beta_u;
         delta_tl = C*delta_x + rho_l;
         delta_tu = -C*delta_x + rho_u;
         delta_wl = -Sl1*Wl*delta_sl -wl;
         delta_wu = -Su1*Wu*delta_su - wu;
         delta_zl = -Tl1*Zl*delta_tl - zl;
         delta_zu = -Tu1*Zu*delta_tu - zu;

        % Determine primal step length alpha_pri_aff
            step = [delta_sl; delta_su; delta_tl; delta_tu];
            curr = [sl; su; tl; tu];
            index = find(step < 0);
            if isempty(index)
                alpha_pri = 1;
            else
                alpha_pri = eta*min(curr(index)./(-step(index)));
                if alpha_pri > 1
                    alpha_pri = 1;
                end
            end
        
        % Determine dual step length alpha_dual_aff
            step = [delta_zl; delta_zu; delta_wl; delta_wu];
            curr = [zl; zu; wl; wu];
            index = find(step < 0);
            if isempty(index)
                alpha_dual = 1;
            else
                alpha_dual = eta*min(curr(index)./(-step(index)));
                if alpha_dual > 1
                    alpha_dual = 1;
                end
            end

       eta = eta + (1-eta)/2;

        % Calculate next iterate
        x = x + alpha_pri * delta_x;
        sl = sl + alpha_pri * delta_sl;
        su = su + alpha_pri * delta_su;
        tl = tl + alpha_pri * delta_tl;
        tu = tu + alpha_pri * delta_tu;
        y = y + alpha_dual * delta_y;
        wl = wl + alpha_dual * delta_wl;
        wu = wu + alpha_dual * delta_wu;
        zl = zl + alpha_dual * delta_zl;
        zu = zu + alpha_dual * delta_zu;
        
        % Convergence
        if mu < 10^(-15)
            break
        end
    end
end