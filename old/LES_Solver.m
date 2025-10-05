function [step] = LES_Solver(dim, obj, equ, ineq, help, curr_it)
%   Solve the required linear equation system 
%   for Interior Point
    
    mat = [help.Hphi         equ.A'           ineq.C';
           equ.A     zeros(dim.m)      zeros(m,r);
           ineq.C  zeros(dim.r, dim.m)   -help.PSI1   ];

    nu = - obj.c + equ.A'*curr_it.y + ineq.C'*(curr_it.zl - curr_it.zu) ...
    - help.bound_xl*help.Sl1*help.Wl*help.betal + help.bound_xu*help.Su1*help.Wu*help.betau;

    kappa = equ.A*curr_it.x -equ.b;

    xi = -curr_it.zl - help.bound_cl*help.Tl1*help.Zl*help.rhol ...
    + curr_it.zu + help.bound_cu*help.Tu1*help.Zu*help.rhou;
    
    omega_3 = help.psi1*xi;
    omega = [nu; -kappa; omega_3];

    % calculate affine step
    sol = mat\omega;
    step = struct();
    step.del_x = sol(1:n);
    
    step.del_sl = help.bound_xl*(step.del_x + help.betal);
    step.del_su = help.bound_xu*(-step.del_x + help.betau);
    step.del_wl = help.bound_xl*(-curr_it.wl -help.Sl1*help.Wl*step.del_sl);
    step.del_wu = help.bound_xu*(-curr_it.wu - help.Su1*help.Wu*step.del_su);
        
    step.del_tl = help.bound_cl*(ineq.C*step.del_x + help.rhol);
    step.del_tu = help.bound_cu*(-ineq.C*step.del_x + help.rhou);
    step.del_zu = help.bound_cu*(-curr_it.zu -help.Tu1*help.Zu*step.del_tu);
    step.del_zl = help.bound_cl*(-curr_it.zl - help.Tl1*help.Zl*step.del_tl);

end
