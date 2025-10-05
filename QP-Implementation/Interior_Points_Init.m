function [it] = Interior_Points_Init(nlp,dim, rnd)
% Calculate initial iterate of Interior Point QP Method 
% return random value for x if desired
% struct it contains x, sl,su,tl,tu,y,wl,wu,zl,zu and
% bound matrices bound_xl, bound_xu, bound_cl, bound_cu

    % problem dimension
    n = dim.n;
    m = dim.m;
    r = dim.r;

    %helping variables
    en = ones(n,1);
    er = ones(r,1);
    
    % set Lagrange multipliers wl, wu 
    it.wl = ones(n,1);
    it.wu = ones(n,1);

    % activate cutest
    p = cutest_setup;

    % set Lagrange multiplier zl to recommended value
    if r == 0
        it.zl = ones(r,1);
    else
        it.zl = max(0.1*ones(r,1), p.v(m+1:m+r));
    end
    it.zu = ones(r,1);
    
    % Set variable x to recommended value or random value if indicated
    it.x = p.x;
    it.y = p.v(1:m,1);
    cutest_terminate
    if rnd
        it.zl = ones(r,1);
        it.x = rand(n,1);
        it.y = rand(m,1);
    end

    % cancel lines where x-component has no upper or no lower bound
    % return bound matrices in sparse format
    index_l = ones(1,n);
    index_u = ones(1,n);
    for i = 1:n
        if nlp.xl(i) < -10^7
            nlp.xl(i) = 0;
            it.wl(i) = 0; 
            index_l(i) = 0;
        end
        if nlp.xu(i) > 10^7
            nlp.xu(i) = 0;
            it.wu(i) = 0;
            index_u(i) = 0;
        end
    end
    index_l = find(index_l);
    index_u = find(index_u);
    it.bound_xl = sparse(index_l, index_l, ones(size(index_l)), n, n);
    it.bound_xu = sparse(index_u, index_u, ones(size(index_u)), n, n);

    it.sl = it.bound_xl*ones(n,1);
    it.su = it.bound_xu*ones(n,1);

    index_cl = ones(1,r);
    index_cu = ones(1,r);
    for i = 1:r
        if nlp.cl(i) < -10^3
            nlp.cl(i) = 0;
            it.zl(i) = 0;
            index_cl(i) = 0;

        end
        if nlp.cu(i) > 10^3
            nlp.cu(i) = 0;
            it.zu(i) = 0;
            index_cu(i) = 0;
        end
    end
    index_cl = find(index_cl);
    index_cu = find(index_cu);
    it.bound_cl = sparse(index_cl, index_cl, ones(size(index_cl)), r, r);
    it.bound_cu = sparse(index_cu, index_cu, ones(size(index_cu)), r, r);

    it.tl = it.bound_cl*ones(r,1);
    it.tu = it.bound_cu*ones(r,1);

    %prepare equation system
    help = helpers(dim, nlp, it);

    % helping variables
    Hphi = nlp.H + help.PHI;
    xi = -it.bound_cl*(it.zl + help.Tl1*help.Zl*help.rho_l) + it.bound_cu*(it.zu + help.Tu1*help.Zu*help.rho_u);
    foo = help.PSI1*xi;

    % calculate affine step with reduced system
    mat = [Hphi nlp.A' nlp.C';
        nlp.A sparse(m,m + r);
        nlp.C sparse(r,m) help.PSI1];
    omega = [-nlp.H*it.x - nlp.c + nlp.A'*it.y + ...
        nlp.C'*it.bound_cl*it.zl - nlp.C'*it.bound_cu*it.zu ...
        - it.bound_xl*help.Sl1*help.Wl*help.beta_l + it.bound_xu*help.Su1*help.Wu*help.beta_u;
        -(nlp.A*it.x-nlp.b);
        foo];
    sol = mat\omega;

    % If mat is singular, compute minimum value of 
    % |mat*del - omega|
    if any(isnan(sol))
        sol = lsqminnorm(mat, omega);
    end
    
    % affine step
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
    
    % duality measure
    it.mu = (it.sl'*it.wl + it.su'*it.wu + it.tl'*it.zl + it.tu'*it.zu)/(2*(n+r));
end
