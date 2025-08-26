function [it,nlp,dim] = Interior_gen_Init()
    
    p = cutest_setup;
    dim.n = p.n;
    dim.m = p.m;
    dim.r = 0;
    
    % compute correct dimensions for m and r
    for i = 1:p.m
        if p.cl(i) || p.cu(i)
        dim.r = dim.r + 1;
        dim.m = dim.m - 1;
        end
    end
    
    n = dim.n;
    m = dim.m;
    r = dim.r;
    en = ones(n,1);
    er = ones(r,1);
    
    %Determine which lines of the constraints are inequalities, which
    %equalities
    j = 0;
    k = 0;
    nlp.equ = zeros(dim.m,1);
    nlp.inequ = zeros(dim.r,1);
    for i = 1:p.m
        if p.cl(i) || p.cu(i)
            k = k+1;
            nlp.inequ(k) = i;
        else
            j = j+1;
            nlp.equ(j) = i;
        end
    end

    % determine all general boundaries for x
    nlp.xl = p.bl;
    nlp.xu = p.bu;
    it.bound_xl = eye(dim.n);
    it.bound_xu = eye(dim.n);

    for i = 1:dim.n
        if nlp.xl(i) < -10^3
            it.bound_xl(i,i) = 0;
            nlp.xl(i) = 0;
        end
        if nlp.xu(i) > 10^3
            it.bound_xu(i,i) = 0;
            nlp.xu(i) = 0;
        end
    end

    %determine the boundaries for Cx
    it.bound_cl = eye(dim.r);
    it.bound_cu = eye(dim.r);
    for i = 1:dim.r
        if p.cl(i) < -10^3
            it.bound_cl(i,i) = 0;
        end
        if p.cu(i) > 10^3
            it.bound_cu = 0;
        end
    end
       
    it.x = ones(dim.n,1);
    it.y = ones(dim.m,1);
    it.sl = it.bound_xl*ones(dim.n,1);
    it.su = it.bound_xu*ones(dim.n,1);
    it.wl = it.bound_xl*ones(dim.n,1);
    it.wu = it.bound_xu*ones(dim.n,1);
    it.tl = it.bound_cl*ones(dim.r,1);
    it.tu = it.bound_cu*ones(dim.r,1);
    it.zl = it.bound_cl*ones(dim.r,1);
    it.zu = it.bound_cu*ones(dim.r,1);

    [nlp] = cutest_iterate(it, nlp, dim,p);
    [help] = helpers_nlp(dim, nlp, it,p);

    mat = [nlp.H + help.PHI   nlp.A'  nlp.C';
           nlp.A   zeros(dim.m)    zeros(dim.m,r);
           nlp.C  zeros(r,dim.m)   -help.PSI1];

    xi = -it.bound_cl*(it.zl + help.Tl1*help.Zl*help.rho_l) + it.bound_cu*(it.zu + help.Tu1*help.Zu*help.rho_u);
    kappa = -nlp.c_e;
    nu = - nlp.grad + nlp.A'*it.y + nlp.C'*(it.zl - it.zu) -it.bound_xl*help.Sl1*help.Wl*help.beta_l ...
            + it.bound_xu*help.Su1*help.Wu*help.beta_u;
    sol = mat\[nu;kappa;xi];

        del_x = sol(1:n);
        del_z = -sol(n+m+1:n+m+r);

        % calculate all variables of expanded system
        del_sl = it.bound_xl*(del_x + help.beta_l);
        del_su = it.bound_xu*(-del_x + help.beta_u);
        del_wl = it.bound_xl*(-it.wl -help.Sl1*help.Wl*del_sl);
        del_wu = it.bound_xu*(-it.wu - help.Su1*help.Wu*del_su);

        del_tl = it.bound_cl*(nlp.C*del_x + help.rho_l);
        del_tu = it.bound_cu*(-nlp.C*del_x + help.rho_u);
        del_zl = -it.bound_cl*(help.Tl1*help.Zl*del_tl + it.zl);
        del_zu = it.bound_cu*(-del_z + del_zl);


        % ensure positivity of s, w, t and z (interiority)
        it.sl = it.bound_xl*max(abs(it.sl + del_sl),en);
        it.su = it.bound_xu*max(abs(it.su + del_su),en);

        it.wl = it.bound_xl*max(abs(it.wl + del_wl),en);
        it.wu = it.bound_xu*max(abs(it.wu + del_wu),en);

        it.tl = it.bound_cl*max(abs(it.tl + del_tl),er);
        it.tu = it.bound_cu*max(abs(it.tu + del_tu),er);

        it.zl = it.bound_cl*max(abs(it.zl + del_zl),er);
        it.zu = it.bound_cu*max(abs(it.zu + del_zu),er);


    cutest_terminate;
end