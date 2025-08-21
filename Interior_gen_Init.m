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
    
    %Determine which lines of the constraints are inequalities, which
    %equalities
    j = 0;
    k = 0;
    nlp.equ = zeros(dim.m,1);
    nlp.inequ = zeros(dim.r,1);
    nlp.cl = zeros(dim.r,1);
    nlp.cu = zeros(dim.r,1);
    for i = 1:p.m
        if any([p.cl(i) p.cu(i)])
            k = k+1;
            nlp.inequ(k) = i;
            nlp.cl(k) = p.cl(i);
            nlp.cu(k) = p.cu(i);
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
        if nlp.cl(i) < -10^3
            it.bound_cl(i,i) = 0;
        end
        if nlp.cu(i) > 10^3
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

    cutest_terminate;
end