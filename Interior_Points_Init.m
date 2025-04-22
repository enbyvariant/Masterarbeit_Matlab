function [x, sl, su, y, wl, wu, bound] = Interior_Points_Init(H, A, b, c, xl, xu, x, y)
% Calculate initial iterate of Interior Point QP Method 

    n = size(A, 2);
    m = size(A, 1);
    e = ones(n,1);

    sl = NaN;
    su = NaN;
    wl = NaN;
    wu = NaN;

    % check if x has upper or lower bound
    bound = "";
    if xl > -10^3
        bound = bound + "l";
        sl = rand(n,1);
        wl = rand(n,1);
        Sl1 = diag(ones(n,1)./sl);
        Wl = diag(wl);
    end
    if xu < 10^3
        bound = bound + "u";
        su = rand(n,1);
        wu = rand(n,1);
        Su1 = diag(ones(n,1)./su);
        Wu = diag(wu);
    end

    % calculate affine steps according to which bound of x exists
    switch bound
        case "lu"
            mat = [H zeros(n,2*n) A' eye(n) -eye(n);
                zeros(n) Sl1*Wl zeros(n,n+m) -eye(n) zeros(n); 
                zeros(n,2*n) Su1*Wu zeros(n, m+n) -eye(n); 
                A zeros(m,m+4*n);
                eye(n) -eye(n) zeros(m+3*n); 
                -eye(n) zeros(n) -eye(n) zeros(m+2*n)];
            omega = [H*x + c - A'*y -wl + wu; Wl*e; Wu*e; A*x - b; (x -xl -sl); (-x +xu -su)];
            sol = linsolve(mat, -omega);
            del_x = sol(1:n);
            del_sl = sol(n+1:2*n);
            del_su = sol(2*n+1:3*n);
            del_wl = sol(3*n+m+1:4*n+m);
            del_wu = sol(4*n+m+1:5*n+m);

            sl = max(abs(sl + del_sl),e);
            su = max(abs(su + del_su),e);
            wl = max(abs(wl + del_wl),e);
            wu = max(abs(wu + del_wu),e);
        case "l"
            mat = [H zeros(n,n) A' eye(n);
                zeros(n) Sl1*Wl zeros(n,m) -eye(n); 
                A zeros(m,m+2*n);
                eye(n) -eye(n) zeros(m+n)];
            omega = [H*x + c - A'*y -wl; Wl*e; A*x -b; x-xl-sl];
            sol = linsolve(mat, -omega);
            del_x = sol(1:n);
            del_sl = sol(n+1:2*n);
            del_wl = sol(2*n+m+1:3*n+m);

            sl = max(abs(sl + del_sl),e);
            wl = max(abs(wl + del_wl),e);
        case "u"
            mat = [H zeros(n,n) A' -eye(n);
                zeros(n) Su1*Wu zeros(n,m) -eye(n); 
                A zeros(m,m+2*n);
                -eye(n) eye(n) zeros(m+n)];
            omega = [H*x + c - A'*y -wu; Wu*e; A*x -b; -x+xu-su];
            sol = linsolve(mat, -omega);
            del_x = sol(1:n);
            del_su = sol(n+1:2*n);
            del_wu = sol(2*n+m+1:3*n+m);

            su = max(abs(su + del_su),e);
            wu = max(abs(wu + del_wu),e);
        otherwise
            mat = [H A'; A zeros(m)];
            omega = [H*x + c - A'*y; A*x - b];
            sol = linsolve(mat, -omega);
            del_x = sol(1:n);
    end
    % ensure positivity of x, s and w (interiority)
    x = max(abs(x + del_x),e);

end