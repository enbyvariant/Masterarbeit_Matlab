function [help] = helpers_nlp(dim, nlp, it,~)
%   compute helping variables
    
    
    help.beta_l = it.bound_xl*(it.x - nlp.xl - it.sl);
    help.beta_u = it.bound_xu*(-it.x + nlp.xu - it.su);
    help.rho_l = it.bound_cl*(nlp.c_i - it.tl);
    help.rho_u = it.bound_cu*(-nlp.c_i - it.tu);

    gen_i = 1:dim.n;
    gen_j = 1:dim.n;
    Slv = zeros(dim.n,1);
    Suv = zeros(dim.n,1);    

    for i = 1:dim.n
        if it.sl(i) ~= 0
            Slv(i) = 1/it.sl(i);
        else
            Slv(i) = 0;
        end
        if it.su(i) ~= 0
            Suv(i) = 1/it.su(i);
        else
            Suv(i) = 0;
        end
    end
    help.Sl1 = sparse(gen_i, gen_j, Slv, dim.n, dim.n);
    help.Su1 = sparse(gen_i, gen_j, Suv, dim.n, dim.n);

    help.Wl = sparse(gen_i, gen_j, it.wl);
    help.Wu = sparse(gen_i, gen_j, it.wu);
            
    gen_i = 1:dim.r;
    gen_j = 1:dim.r;
    Tlv = zeros(dim.r,1);
    Tuv = zeros(dim.r,1);    
    for i = 1:dim.r
        if it.tl(i) ~= 0
            Tlv(i) = 1/it.tl(i);
        else
            Tlv(i) = 0;
        end
        if it.tu(i) ~= 0
            Tuv(i) = 1/it.tu(i);
        end
    end
    help.Tl1 = sparse(gen_i,gen_j, Tlv, dim.r, dim.r);
    help.Tu1 = sparse(gen_i, gen_j, Tuv, dim.r,dim.r);

    help.Zl = sparse(gen_i, gen_j, it.zl, dim.r, dim.r);
    help.Zu = sparse(gen_i, gen_j, it.zu, dim.r, dim.r);     

    help.PHI = sparse(it.bound_xl*help.Sl1*help.Wl + it.bound_xu*help.Su1*help.Wu);
    PSI = it.bound_cl*help.Tl1*help.Zl + it.bound_cu*help.Tu1*help.Zu;
    Pv = zeros(dim.r,1);
        for i = 1:dim.r
            if PSI(i,i) ~= 0
                Pv(i) = 1/PSI(i,i);
            else
                Pv(i) = 0;
            end
        end
    help.PSI1 = sparse(gen_i, gen_j, Pv, dim.r, dim.r);

end