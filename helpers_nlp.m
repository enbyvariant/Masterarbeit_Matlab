function [help] = helpers_nlp(dim, nlp, it,p)
%   compute helping variables
    
    
    help.beta_l = it.x - nlp.xl - it.sl;
    help.beta_u = -it.x + nlp.xu - it.su;
    help.rho_l = (nlp.c_i - it.tl);
    help.rho_u = -nlp.c_i - it.tu;

    help.Sl1 = sparse(zeros(dim.n));
    help.Su1 = sparse(zeros(dim.n));
    for i = 1:dim.n
        if it.sl(i) ~= 0
            help.Sl1(i,i) = 1/it.sl(i);
        end
        if it.su(i) ~= 0
            help.Su1(i,i) = 1/it.su(i);
        end
    end
    help.Wl = sparse(diag(it.wl));
    help.Wu = sparse(diag(it.wu));
            
    help.Tl1 = sparse(zeros(dim.r));
    help.Tu1 = sparse(zeros(dim.r));
    for i = 1:dim.r
        if it.tl(i) ~= 0
            help.Tl1(i,i) = 1/it.tl(i);
        end
        if it.tu(i) ~= 0
            help.Tu1(i,i) = 1/it.tu(i);
        end
    end
    help.Zl = sparse(diag(it.zl));
    help.Zu = sparse(diag(it.zu));               
    help.PHI = it.bound_xl*help.Sl1*help.Wl + it.bound_xu*help.Su1*help.Wu;
    PSI = it.bound_cl*help.Tl1*help.Zl + it.bound_cu*help.Tu1*help.Zu;
    help.PSI1 = zeros(dim.r);
        for i = 1:dim.r
            if PSI(i,i) ~= 0
            help.PSI1(i,i) = 1/PSI(i,i);
            end
        end

end