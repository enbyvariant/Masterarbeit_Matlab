function [x, sl, su, y, wl, wu, bound_l, bound_u] = Interior_Points_Init(H, A, b, c, xl, xu, x, y)
% Calculate initial iterate of Interior Point QP Method 

    n = size(A, 2);
    m = size(A, 1);
    e = ones(n,1);

    sl = zeros(n,1);
    su = zeros(n,1);
    wl = zeros(n,1);
    wu = zeros(n,1);

    % check if each x-component has upper or lower bound
    bound_l = zeros(n);
    bound_u = zeros(n);
    for i = 1:n
        if xl(i) > -10^3
            bound_l(i,i) = 1;
            sl(i) = rand(1);
            wl(i) = rand(1);
        end
        if xu < 10^3
            bound_u(i,i) = 1;
            su(i) = rand(1);
            wu(i) = rand(1);
        end
    end
    Sl1 = diag(ones(n,1)./sl);
    Wl = diag(wl);
    Su1 = diag(ones(n,1)./su);
    Wu = diag(wu);

    % calculate affine step according to which bound of x exists

    mat = [H zeros(n,2*n) A' bound_l*eye(n) -bound_u*eye(n);
        zeros(n) bound_l*Sl1*Wl zeros(n,n+m) -bound_l*eye(n) zeros(n);
        zeros(n,2*n) bound_u*Su1*Wu zeros(n, m+n) -bound_u*eye(n);
        A zeros(m,m+4*n);
        bound_l*eye(n) -bound_l*eye(n) zeros(m+3*n);
        -bound_u*eye(n) zeros(n) -bound_u*eye(n) zeros(m+2*n)];
    omega = [H*x + c - A'*y - wl + wu; Wl*e; Wu*e; A*x - b; bound_l*(x -xl -sl); bound_u*(-x +xu -su)];
    sol = omega\mat;
    
    del_x = sol(1:n);
    del_sl = sol(n+1:2*n);
    del_su = sol(2*n+1:3*n);
    del_wl = sol(3*n+m+1:4*n+m);
    del_wu = sol(4*n+m+1:5*n+m);
    
    % ensure positivity of x, s and w (interiority)
    sl = max(abs(sl + del_sl),e);
    su = max(abs(su + del_su),e);

    wl = max(abs(wl + del_wl),e);
    wu = max(abs(wu + del_wu),e);

    x = max(abs(x + del_x),e);
end