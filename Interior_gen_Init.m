function [it,nlp,dim, f] = Interior_gen_Init(cue,p)
% Calculate initial iterate of Interior Point NLP Method
% Calculate initial problem data nlp
% struct it contains x, sl,su,tl,tu,y,wl,wu,zl,zu and
% bound matrices bound_xl, bound_xu, bound_cl, bound_cu
% struct nlp contains H, obj, xl, xu, cl, cu, A, C
   
    dim.n = p.n;
    dim.m = p.m;
    dim.r = 0;
    f = 0;
    % uncomment if only smaller problems are desired
    % if dim.n > 20
    %     f = 1;
    %     nlp = [];
    %     it = [];
    %     fprintf(1, 'n = %i\n', dim.n);
    %     fprintf(1, 'm = %i\n', dim.m);
    %     return
    % end

    % compute correct dimensions for m and r
    for i = 1:p.m
        if p.cl(i) || p.cu(i)
        dim.r = dim.r + 1;
        dim.m = dim.m - 1;
        end
    end
    
    % helping variables
    n = dim.n;
    m = dim.m;
    r = dim.r;
    en = ones(n,1);
    er = ones(r,1);
    
    % test for missing constraints
    free = 1;
    for i=1:n
       if p.bl(i) > -10^3 || p.bu(i) < 10^3
           free = 0;
           continue
       end
    end
    if free 
       p.bl = -10*ones(n,1);
       p.bu = 10*ones(n,1);
    end
    
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
    gen_i = 1:dim.n;
    gen_j = 1:dim.n;
    index_xl = zeros(dim.n,1);
    index_xu = zeros(dim.n,1);

    for i = 1:dim.n
        if nlp.xl(i) > -10^7
            index_xl(i) = 1;
        end
        if nlp.xu(i) < 10^7
            index_xu(i) = 1;
        end
    end
    it.bound_xl = sparse(gen_i, gen_j, index_xl, dim.n, dim.n);
    it.bound_xu = sparse(gen_i, gen_j, index_xl, dim.n, dim.n);

    %determine the boundaries for Cx
    index_cl = ones(dim.r,1);
    index_cu = ones(dim.r,1);
    for i = 1:dim.r
        if p.cl(i) < -10^7
            index_cl(i) = 0;
        end
        if p.cu(i) > 10^7
            index_cu(i) = 0;
         end
    end
    gen_i = 1:dim.r;
    gen_j = 1:dim.r;
    it.bound_cl = sparse(gen_i, gen_j, index_cl, dim.r,dim.r);
    it.bound_cu = sparse(gen_i, gen_j, index_cu, dim.r,dim.r);
       
    it.x = sparse(p.x);
    it.y = sparse(p.v(1:dim.m));
    if dim.m == 0
        it.y = zeros(dim.m,1);
    end
    it.sl = sparse(it.bound_xl*ones(dim.n,1));
    it.su = sparse(it.bound_xu*ones(dim.n,1));
    it.wl = sparse(it.bound_xl*ones(dim.n,1));
    it.wu = sparse(it.bound_xu*ones(dim.n,1));
    it.tl = sparse(it.bound_cl*ones(dim.r,1));
    it.tu = sparse(it.bound_cu*ones(dim.r,1));
    if dim.r ~= 0
        it.zl = sparse(it.bound_cl*p.v(dim.m+1:dim.m + dim.r));
        it.zu = sparse(it.bound_cu*p.v(dim.m+1:dim.m + dim.r));
    else
        it.zl = zeros(dim.r,1);
        it.zu = zeros(dim.r,1);
    end

    it.mu = (it.sl'*it.wl + it.su'*it.wu + it.tl'*it.zl + it.tu'*it.zu)/(2*dim.r+2*dim.n);
    
    [nlp] = cutest_iterate(it, nlp, dim,p);
    [help] = helpers_nlp(dim, nlp, it,p);
    if cue
        [gamma] = matrix_factors(dim, nlp);
    else
        gamma = 0;
    end
    Hphi = nlp.H + help.PHI + gamma*eye(n);
    mat = [Hphi   nlp.A'  nlp.C';
           nlp.A   sparse(m,m)    sparse(dim.m,r);
           nlp.C  sparse(r,dim.m)   -help.PSI1];

    xi = -it.bound_cl*(it.zl + help.Tl1*help.Zl*help.rho_l) + it.bound_cu*(it.zu + help.Tu1*help.Zu*help.rho_u);
    kappa = -nlp.c_e;
    nu = - nlp.grad + nlp.A'*it.y + nlp.C'*(it.zl - it.zu) -it.bound_xl*help.Sl1*help.Wl*help.beta_l ...
            + it.bound_xu*help.Su1*help.Wu*help.beta_u;
    sol = mat\[nu;kappa;help.PSI1*xi];
    if any(isnan(sol))
        sol = lsqminnorm(mat, [nu;kappa;help.PSI1*xi]);
    end

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
        it.sl = it.bound_xl*max(abs(it.sl + del_sl),0.1*en);
        it.su = it.bound_xu*max(abs(it.su + del_su),0.1*en);

        it.wl = it.bound_xl*max(abs(it.wl + del_wl),0.1*en);
        it.wu = it.bound_xu*max(abs(it.wu + del_wu),0.1*en);

        it.tl = it.bound_cl*max(abs(it.tl + del_tl),0.1*er);
        it.tu = it.bound_cu*max(abs(it.tu + del_tu),0.1*er);

        it.zl = it.bound_cl*max(abs(it.zl + del_zl),0.1*er);
        it.zu = it.bound_cu*max(abs(it.zu + del_zu),0.1*er);

end