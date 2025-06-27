function [x,y,sl,su,wl,wu,tl, tu,zl,zu,obj] = Interior_Points_NLP(max_iter)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    
    [x,y,sl,su,wl,wu,tl,tu,zl,zu,bound_xl,bound_xu,bound_cl,bound_cu, equ, inequ,n, m, r] = Interior_gen_Init();
    iter = 0;
    while 1
        if iter > max_iter
            break
        else
            iter = iter + 1;
        end

        % calculate variables for new iterate
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
            PSI1(i,i) =
        end
        Zl = diag(zl);
        Zu = diag(zu);

        PHI = Sl1*Wl+Su1*Wu;
        PSI = Tl1*Zl+Tu1*Zu;
        PSI1 = zeros(r);
        for i = 1:r
            PSI1(i,i) = 1/PSI(i,i);
        end
        [grad, H, cons_e, cons_i, cons_e_grad, cons_i_grad, obj] = cutest_iterate(x,y,zl,zu,equ,inequ);

        mat = [H + PHI   cons_e_grad'  cons_i_grad';
             cons_e_grad   zeros(m)    zeros(m,r);
             cons_i_grad  zeros(r,m)   -PSI1];

    end

end