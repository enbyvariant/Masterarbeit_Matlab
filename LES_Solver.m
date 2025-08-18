function [step] = LES_Solver(dim, nlp, help, it, sigma, step_a)
%   Solve the required linear equation system 
%   for Interior Point

    en = ones(dim.n,1);
    er = ones(dim.r,1);
    %disp(sigma);
    phi_l = it.bound_xl*(it.wl - help.Sl1*(sigma*it.mu*en - diag(step_a.del_sl)*step_a.del_wl));
    phi_u = it.bound_xu*(it.wu - help.Su1*(sigma*it.mu*en - diag(step_a.del_su)*step_a.del_wu));
    psi_l = it.bound_cl*(it.zl - help.Tl1*(sigma*it.mu*er - diag(step_a.del_tl)*step_a.del_zl));
    psi_u = it.bound_cu*(it.zu - help.Tu1*(sigma*it.mu*er - diag(step_a.del_tu)*step_a.del_zu));
    
    % Set up equation system
    mat = [help.Hphi       nlp.A'        nlp.C';
           nlp.A     zeros(dim.m)      zeros(dim.m,dim.r);
           nlp.C  zeros(dim.r, dim.m)   -help.PSI1   ];
    nabla_L = nlp.c - nlp.A'*it.y - nlp.C'*(it.bound_cl*it.zl - ...
        it.bound_cu*it.zu) - it.bound_xl*it.wl + it.bound_xu*it.wu;
    omega_1 = -(nabla_L + it.bound_xl*(phi_l + help.Sl1*help.Wl*help.beta_l) - ...
    it.bound_xu*(phi_u + help.Su1*help.Wu*help.beta_u));
    xi = -psi_l + psi_u - it.bound_cl*help.Tl1*help.Zl*help.rho_l + it.bound_cu*help.Tu1*help.Zu*help.rho_u;
    omega_3 = help.PSI1*xi;
    omega = [omega_1;-help.kappa;omega_3];
    
    % nu = nabla_L + phi_l + it.bound_xl*help.Sl1*help.Wl*help.beta_l - phi_u - it.bound_xu*help.Su1*help.Wu*help.beta_u;
    % 
    % xi = psi_l + it.bound_cl*help.Tl1*help.Zl*help.rho_l ...
    % - psi_u - it.bound_cu*help.Tu1*help.Zu*help.rho_u;
    % 
    % omega_3 = help.PSI1*xi;
    % omega = [-nu; -help.kappa; -omega_3];

    % calculate affine step
    sol = mat\omega;
    step = struct();
    step.del_x = sol(1:dim.n);
    disp(step.del_x)
    step.del_y = -sol(dim.n+1:dim.n+dim.m);
    
    step.del_sl = it.bound_xl*(step.del_x + help.beta_l);
    step.del_su = it.bound_xu*(-step.del_x + help.beta_u);
    step.del_wl = it.bound_xl*(-phi_l -help.Sl1*help.Wl*step.del_sl);
    step.del_wu = it.bound_xu*(-phi_u - help.Su1*help.Wu*step.del_su);
        
    step.del_tl = it.bound_cl*(nlp.C*step.del_x + help.rho_l);
    step.del_tu = it.bound_cu*(-nlp.C*step.del_x + help.rho_u);
    step.del_zu = it.bound_cu*(-psi_l -help.Tu1*help.Zu*step.del_tu);
    step.del_zl = it.bound_cl*(-psi_u - help.Tl1*help.Zl*step.del_tl);

end