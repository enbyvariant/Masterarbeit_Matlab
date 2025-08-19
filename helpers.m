function [help] = helpers(dim, nlp, it)
%   UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    help.beta_l = it.x - nlp.xl - it.sl;
    help.beta_u = -it.x + nlp.xu - it.su;
    help.rho_l = nlp.C*x - nlp.cl - it.tl;
    help.rho_u = -nlp.C*x + nlp.cu - it.tu;

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
end