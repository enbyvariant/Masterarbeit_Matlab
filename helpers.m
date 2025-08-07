function [help] = helpers(dim, obj, ineq, equ, x_b, it)
%   UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    help.beta_l = it.x - x_b.xl - it.sl;
    help.beta_u = -it.x + x_b.xu - it.su;
    help.rho_l = ineq.C*x - ineq.cl - it.tl;
    help.rho_u = -ineq.C*x + ineq.cu - it.tu;

    Sl1 = zeros(n);
    Su1 = zeros(n);
    for i = 1:n
        if sl(i) ~= 0
            Sl1(i,i) = 1/sl(i);
        end
        if su(i) ~= 0
            Su1(i,i) = 1/su(i);
        end
    end
    Wl = diag(wl);
    Wu = diag(wu);
            
    Tl1 = zeros(r);
    Tu1 = zeros(r);
    for i = 1:r
        if tl(i) ~= 0
            Tl1(i,i) = 1/tl(i);
        end
        if tu(i) ~= 0
            Tu1(i,i) = 1/tu(i);
        end
    end
    Zl = diag(zl);
    Zu = diag(zu);   
end