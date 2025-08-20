function [help] = helpers(dim, nlp, it)
%   UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    help.beta_l = it.x - nlp.xl - it.sl;
    help.beta_u = -it.x + nlp.xu - it.su;
    help.rho_l = nlp.C*it.x - nlp.cl - it.tl;
    help.rho_u = -nlp.C*it.x + nlp.cu - it.tu;

    help.Sl1 = zeros(dim.n);
    help.Su1 = zeros(dim.n);
    for i = 1:dim.n
        if it.sl(i) ~= 0
            help.Sl1(i,i) = 1/it.sl(i);
        end
        if it.su(i) ~= 0
            help.Su1(i,i) = 1/it.su(i);
        end
    end
    help.Wl = diag(it.wl);
    help.Wu = diag(it.wu);
            
    help.Tl1 = zeros(dim.r);
    help.Tu1 = zeros(dim.r);
    for i = 1:dim.r
        if it.tl(i) ~= 0
            help.Tl1(i,i) = 1/it.tl(i);
        end
        if it.tu(i) ~= 0
            help.Tu1(i,i) = 1/it.tu(i);
        end
    end
    help.Zl = diag(it.zl);
    help.Zu = diag(it.zu);               
    help.PHI = it.bound_xl*help.Sl1*help.Wl + it.bound_xu*help.Su1*help.Wu;
    PSI = it.bound_cl*help.Tl1*help.Zl + it.bound_cu*help.Tu1*help.Zu;
    help.PSI1 = zeros(dim.r);
        for i = 1:dim.r
            if PSI(i,i) ~= 0
            help.PSI1(i,i) = 1/PSI(i,i);
            end
        end

end