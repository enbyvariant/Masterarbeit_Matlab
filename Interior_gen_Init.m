function [x,y,sl,su,wl,wu,tl,tu,zl,zu,bound_xl,bound_xu,bound_cl,bound_cu, equ, inequ,n, m, r] = Interior_gen_Init()
    
    p = cutest_setup;
    n = p.n;
    m = prob.m;
    r = 0;
    
    % compute correct dimensions for m and r
    for i = 1:m
        if prob.cl(i) || prob.cu(i)
        r = r + 1;
        m = m - 1;
        end
    end

    j = 0;
    k = 0;
    equ = zeros(m,1);
    inequ = zeros(r,1);
    cl = zeros(r,1);
    cu = zeros(r,1);
    for i = 1:p.m
        if p.cl(i) == 0 && p.cu(i) == 0
            j = j+1;
            equ(j) = i;
        else
            k = k+1;
            inequ(k) = i;
            cl(k) = p.cl(i);
            cu(k) = p.cu(i);
        end
    end
    xl = p.bl;
    xu = p.bu;
    bound_xl = eye(n);
    bound_xu = eye(n);

    for i = 1:n
        if xl(i) < -10^3
            bound_xl(i,i) = 0;
            xl(i) = 0;
        end
        if xu(i) > 10^3
            bound_xu(i,i) = 0;
            xu(i) = 0;
        end
    end

    bound_cl = eye(r);
    bound_cu = eye(r);
    for i = 1:r
        if cl(i) < -10^3
            bound_cl(i,i) = 0;
        end
        if cu(i) > 10^3
            bound_cu = 0;
        end
    end
       
    x = ones(n,1);
    y = ones(m,1);
    sl = bound_xl*ones(n,1);
    su = bound_xu*ones(n,1);
    wl = bound_xl*ones(n,1);
    wu = bound_xu*ones(n,1);
    tl = bound_cl*ones(r,1);
    tu = bound_cu*ones(r,1);
    zl = bound_cl*ones(r,1);
    zu = bound_cu*ones(r,1);

end