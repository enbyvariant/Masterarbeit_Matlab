function [it, dim, nlp] = Interior_gen_Init()
    
    p = cutest_setup;
    dim.n = p.n;
    dim.m = p.m;
    dim.r = 0;
    
    % compute correct dimensions for m and r
    for i = 1:dim.m
        if p.cl(i) || p.cu(i)
        dim.r = dim.r + 1;
        dim.m = dim.m - 1;
        end
    end
    

    j = 0;
    k = 0;
    nlp.equ = zeros(dim.m,1);
    nlp.inequ = zeros(dim.r,1);
    nlp.cl = zeros(dim.r,1);
    nlp.cu = zeros(dim.r,1);
    for i = 1:p.m
        if p.cl(i) == 0 && p.cu(i) == 0
            j = j+1;
            nlp.equ(j) = i;
        else
            k = k+1;
            nlp.inequ(k) = i;
            nlp.cl(k) = p.cl(i);
            nlp.cu(k) = p.cu(i);
        end
    end
    nlp.xl = p.bl;
    nlp.xu = p.bu;
    nlp.bound_xl = eye(dim.n);
    nlp.bound_xu = eye(dim.n);

    for i = 1:dim.n
        if nlp.xl(i) < -10^3
            nlp.bound_xl(i,i) = 0;
            nlp.xl(i) = 0;
        end
        if nlp.xu(i) > 10^3
            nlp.bound_xu(i,i) = 0;
            nlp.xu(i) = 0;
        end
    end

    nlp.bound_cl = eye(dim.r);
    nlp.bound_cu = eye(dim.r);
    for i = 1:dim.r
        if nlp.cl(i) < -10^3
            nlp.bound_cl(i,i) = 0;
        end
        if nlp.cu(i) > 10^3
            nlp.bound_cu = 0;
        end
    end
       
    it.x = ones(dim.n,1);
    it.y = ones(dim.m,1);
    it.sl = nlp.bound_xl*ones(dim.n,1);
    it.su = nlp.bound_xu*ones(dim.n,1);
    it.wl = nlp.bound_xl*ones(dim.n,1);
    it.wu = nlp.bound_xu*ones(dim.n,1);
    it.tl = nlp.bound_cl*ones(dim.r,1);
    it.tu = nlp.bound_cu*ones(dim.r,1);
    it.zl = nlp.bound_cl*ones(dim.r,1);
    it.zu = nlp.bound_cu*ones(dim.r,1);

    cutest_terminate;
end