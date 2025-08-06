
% Implementation of the Interior Points Method for LPs

function [x,sl,su,tl,tu,y,wl,wu,zl,zu,mu, iterations, data] = Interior_Points_LP(iter, A, C, cl, cu, xl, xu, b, c)
    
    % Initial values
    n = size(A, 2);
    m = size(A, 1);
    r = size(C, 1);
    eta = .985;
    en = ones(n,1);
    er = ones(r,1);

    % Compute starting point
    [x, sl, su, tl, tu, y, wl, wu, zl, zu, bound_xl, bound_xu, bound_cl, bound_cu] = Init_LP(A, C, cl, cu, xl, xu, b, c);
    data = [];
    p = cutest_setup;

    % Main Loop
    for iterations = 1:iter
       %  if r == 0
       %  % Inequalities do not exist
       % 
       % 
       %  % Update variables
       %  Sl1 = zeros(n);
       %  Su1 = zeros(n);
       %  for i = 1:n
       %      if sl(i) ~= 0
       %          Sl1(i,i) = 1/sl(i);
       %      end
       %      if su(i) ~= 0
       %          Su1(i,i) = 1/su(i);
       %      end
       %  end
       %  Wl = diag(wl);
       %  Wu = diag(wu);
       % 
       %  % Define helping variables
       %  PHI = bound_xl*Sl1*Wl + bound_xu*Su1*Wu;
       %  beta_l = bound_xl*(x - xl - sl);
       %  beta_u = bound_xu*(-x + xu -su);
       %  % set up equation system
       %  nu = -(c - A'*y) -bound_xl*Sl1*Wl*beta_l + bound_xu*Su1*Wu*beta_u;
       %  mat = [PHI A';
       %      A zeros(m,m)];
       % 
       %  % Solve for the affine search direction
       %   sol = mat\[nu;-(A*x - b)];
       % 
       %   delta_x = sol(1:n);
       %   delta_sl = bound_xl*(delta_x + beta_l);
       %   delta_su = bound_xu*(-delta_x + beta_u);
       %   delta_wl = bound_xl*(-Sl1*Wl*delta_sl -wl);
       %   delta_wu = bound_xu*(-Su1*Wu*delta_su - wu);
       % 
       %  % Determine affine primal step length alpha_pri_aff
       %      step = [delta_sl; delta_su];
       %      curr = [sl; su];
       %      index = find(step < 0);
       %      if isempty(index)
       %          alpha_pri = 1;
       %      else
       %          alpha_pri = min(curr(index)./(-step(index)));
       %          if alpha_pri > 1
       %              alpha_pri = 1;
       %          end
       %      end
       % 
       %  % Determine affine dual step length alpha_dual_aff
       %      step = [delta_wl; delta_wu];
       %      curr = [wl; wu];
       %      index = find(step < 0);
       %      if isempty(index)
       %          alpha_dual = 1;
       %      else
       %          alpha_dual = min(curr(index)./(-step(index)));
       %          if alpha_dual > 1
       %              alpha_dual = 1;
       %          end
       %      end
       % 
       %  % Calculate duality measure mu
       %  mu = (sl'*wl + su'*wu)/(2*n);
       % 
       %  % Calculate affine duality measure mu_aff
       %  sl_a = sl + alpha_pri*delta_sl;
       %  su_a = su + alpha_pri*delta_su;
       %  wl_a = wl + alpha_dual*delta_wl;
       %  wu_a = wu + alpha_dual*delta_wu;
       %  mu_aff = (sl_a'*wl_a + su_a'*wu_a)/(2*n);
       % 
       %  % Set the centering parameter
       %  sigma = (mu_aff/mu)^3;
       % 
       %  % Calculate the search direction
       %  phi_l = wl - bound_xl*Sl1*(sigma*mu*en - diag(delta_sl)*delta_wl);
       %  phi_u = wu - bound_xu*Su1*(sigma*mu*en - diag(delta_su)*delta_wu);
       %  nu = -(c - A'*y- wl + wu) - phi_l - bound_xl*Sl1*Wl*beta_l + phi_u + bound_xu*Su1*Wu*beta_u;
       %  sol = mat\[nu;-(A*x - b)];
       % 
       %   delta_x = sol(1:n);
       %   delta_y = -sol(n+1:m+n);
       %   delta_sl = bound_xl*(delta_x + beta_l);
       %   delta_su = bound_xu*(-delta_x + beta_u);
       %   delta_wl = bound_xl*(-Sl1*Wl*delta_sl -wl);
       %   delta_wu = bound_xu*(-Su1*Wu*delta_su - wu);
       % 
       %  % Determine primal step length alpha_pri_aff
       %      step = [delta_sl; delta_su];
       %      curr = [sl; su];
       %      index = find(step < 0);
       %      if isempty(index)
       %          alpha_pri = 1;
       %      else
       %          alpha_pri = eta*min(curr(index)./(-step(index)));
       %          if alpha_pri > 1
       %              alpha_pri = 1;
       %          end
       %      end
       % 
       %  % Determine dual step length alpha_dual_aff
       %      step = [delta_wl; delta_wu];
       %      curr = [wl; wu];
       %      index = find(step < 0);
       %      if isempty(index)
       %          alpha_dual = 1;
       %      else
       %          alpha_dual = eta*min(curr(index)./(-step(index)));
       %          if alpha_dual > 1
       %              alpha_dual = 1;
       %          end
       %      end
       % 
       % %eta = eta + (1-eta)/2;
       % 
       %  % Calculate next iterate
       %  x = x + alpha_pri * delta_x;
       %  sl = bound_xl*(sl + alpha_pri * delta_sl);
       %  su = bound_xu*(su + alpha_pri * delta_su);
       %  y = y + alpha_dual * delta_y;
       %  wl = bound_xl*(wl + alpha_dual * delta_wl);
       %  wu = bound_xu*(wu + alpha_dual * delta_wu);
       % 
       %  else
        % Inequalities exist cl <= Cx <= cu
        
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
        PHI = bound_xl*Sl1*Wl + bound_xu*Su1*Wu;
        PSI = bound_cl*Tl1*Zl + bound_cu*Tu1*Zu;
        PSI1 = zeros(r);
        for i = 1:r
            PSI1(i,i) = 1/PSI(i,i);
        end

        beta_l = bound_xl*(x - xl - sl);
        beta_u = bound_xu*(-x + xu -su);
        rho_l = bound_cl*(C*x - cl - tl);
        rho_u = bound_cu*(-C*x + cu - tu);
        
        % set up equation system
        nabla_L = c - A'*y - C'*(zl-zu) - wl - wu;
        nu = (nabla_L) + wl + bound_xl*Sl1*Wl*beta_l - wu - bound_xu*Su1*Wu*beta_u;
        xi = PSI1*(zl + bound_cl*Tl1*Zl*rho_l - zu - bound_cu*Tu1*Zu*rho_u);
        kappa = (A*x - b);
        mat = [PHI A' C';
            A zeros(m,m) zeros(m,r);
            C zeros(r,m) -PSI1];
        omega = [-nu;-kappa;-xi];
        % Solve for the affine search direction
         sol = mat\omega;

         delta_x = sol(1:n);
         delta_sl = bound_xl*(delta_x + beta_l);
         delta_su = bound_xu*(-delta_x + beta_u);
         delta_tl = bound_cl*(C*delta_x + rho_l);
         delta_tu = bound_cu*(-C*delta_x + rho_u);
         delta_wl = bound_xl*(-Sl1*Wl*delta_sl -wl);
         delta_wu = bound_xu*(-Su1*Wu*delta_su - wu);
         delta_zl = bound_cl*(-Tl1*Zl*delta_tl - zl);
         delta_zu = bound_cu*(-Tu1*Zu*delta_tu - zu);
        
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
        phi_l = wl - bound_xl*Sl1*(sigma*mu*en - diag(delta_sl)*delta_wl);
        phi_u = wu - bound_xu*Su1*(sigma*mu*en - diag(delta_su)*delta_wu);
        psi_l = zl - bound_cl*Tl1*(sigma*mu*er - diag(delta_tl)*delta_zl);
        psi_u = zu - bound_cu*Tu1*(sigma*mu*er - diag(delta_tu)*delta_zu);
        nu = nabla_L + phi_l + Sl1*Wl*beta_l - phi_u - Su1*Wu*beta_u;
        xi = PSI1*(psi_l + Tl1*Zl*rho_l - psi_u - Tu1*Zu*rho_u);
        sol = mat\[-nu;-kappa; -xi];
        
         delta_x = sol(1:n);
         delta_y = -sol(n+1:m+n);
         delta_sl = bound_xl*(delta_x + beta_l);
         delta_su = bound_xu*(-delta_x + beta_u);
         delta_tl = bound_cl*(C*delta_x + rho_l);
         delta_tu = bound_cu*(-C*delta_x + rho_u);
         delta_wl = bound_xl*(-Sl1*Wl*delta_sl -wl);
         delta_wu = bound_xu*(-Su1*Wu*delta_su - wu);
         delta_zl = bound_cl*(-Tl1*Zl*delta_tl - zl);
         delta_zu = bound_cu*(-Tu1*Zu*delta_tu - zu);

        % Determine primal step length alpha_pri
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
        
        % Determine dual step length alpha_dual
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

  %     eta = eta + (1-eta)/2;

        % Calculate next iterate
        x = x + alpha_pri * delta_x;
        sl = bound_xl*(sl + alpha_pri * delta_sl);
        su = bound_xu*(su + alpha_pri * delta_su);
        tl = bound_cl*(tl + alpha_pri * delta_tl);
        tu = bound_cu*(tu + alpha_pri * delta_tu);
        y = y + alpha_dual * delta_y;
        wl = bound_xl*(wl + alpha_dual * delta_wl);
        wu = bound_xu*(wu + alpha_dual * delta_wu);
        zl = bound_cl*(zl + alpha_dual * delta_zl);
        zu = bound_cu*(zu + alpha_dual * delta_zu);
            mu_n = (sl'*wl + su'*wu + tl'*zl + tu'*zu)/(2*(n+r));
            obj = cutest_obj(x);

            data = [data; iterations obj mu_n mu-mu_n];
        
        %end
        % Convergence
        if mu < 10^(-15)
            break;
        end
    end
    cutest_terminate
end