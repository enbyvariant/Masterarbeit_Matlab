function [it] = Interior_Points_Init(nlp,dim)
% Calculate initial iterate of Interior Point QP Method 

    n = dim.n;
    m = dim.m;
    r = dim.r;
    en = ones(n,1);
    er = ones(r,1);

    it.wl = ones(n,1);
    it.wu = ones(n,1);

    it.zl = ones(r,1);
    it.zu = ones(r,1);

    it.x = rand(n,1);
    it.y = rand(m,1);

    % cancel lines where x-component has no upper or no lower bound
    index_l = ones(1,n);
    index_u = ones(1,n);
    for i = 1:n
        if nlp.xl(i) < -10^3
            nlp.xl(i) = 0;
            it.wl(i) = 0; 
            index_l(i) = 0;
        end
        if nlp.xu(i) > 10^3
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

    %prepare equation system
    Sl1 = zeros(n);
    Su1 = zeros(n);
    for i = 1:n
        if it.sl(i) ~= 0
            Sl1(i,i) = 1/it.sl(i);
        end
        if it.su(i) ~= 0
            Su1(i,i) = 1/it.su(i);
        end
    end
    Wl = diag(it.wl);
    Wu = diag(it.wu);

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
        Zl = diag(it.zl);
        Zu = diag(it.zu);
        index_cl = find(index_cl);
        index_cu = find(index_cu);
        it.bound_cl = sparse(index_cl, index_cl, ones(size(index_cl)), r, r);
        it.bound_cu = sparse(index_cu, index_cu, ones(size(index_cu)), r, r);

        it.tl = it.bound_cl*ones(r,1);
        it.tu = it.bound_cu*ones(r,1);
        Tl1 = zeros(r);
        Tu1 = zeros(r);
        for i = 1:r
            if it.tl(i) ~= 0
            Tl1(i,i) = 1/it.tl(i);
            end
            if it.tu(i) ~= 0
                Tu1(i,i) = 1/it.tu(i);
            end
        end

        % helping variables
        rhol = it.bound_cl*(nlp.C*it.x-nlp.cl-it.tl);
        rhou = it.bound_cu*(-nlp.C*it.x + nlp.cu - it.tu);
        betal = it.bound_xl*(it.x-nlp.xl-it.sl);
        betau = it.bound_xu*(-it.x+nlp.xu-it.su);
        Hphi = nlp.H + Sl1*Wl + Su1*Wu;
        psi = Tl1*Zl +Tu1*Zu;
        xi = -it.bound_cl*(it.zl + Tl1*Zl*rhol) + it.bound_cu*(it.zu + Tu1*Zu*rhou);
        foo = psi\xi;

        % calculate affine step with reduced system
        mat = [Hphi nlp.A' nlp.C';
            nlp.A zeros(m,m + r);
            nlp.C zeros(r,m) inv(psi)];
        omega = [-nlp.H*it.x - nlp.c + nlp.A'*it.y + nlp.C'*it.bound_cl*it.zl - nlp.C'*it.bound_cu*it.zu - it.bound_xl*Sl1*Wl*betal + it.bound_xu*Su1*Wu*betau;
            -(nlp.A*it.x-nlp.b);
            foo];
        sol = mat\omega;

        del_x = sol(1:n);
        del_z = -sol(n+m+1:n+m+r);

        % calculate all variables of expanded system
        del_sl = it.bound_xl*(del_x + betal);
        del_su = it.bound_xu*(-del_x + betau);
        del_wl = it.bound_xl*(-it.wl -Sl1*Wl*del_sl);
        del_wu = it.bound_xu*(-it.wu - Su1*Wu*del_su);

        del_tl = it.bound_cl*(nlp.C*del_x + rhol);
        del_tu = it.bound_cu*(-nlp.C*del_x + rhou);
        del_zl = -it.bound_cl*(Tl1*Zl*del_tl + it.zl);
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
        
    
end