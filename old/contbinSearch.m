function alpha = contbinSearch(a,b,x,del,tol)
% Compute largest alpha in (a,b], such that 
%   x + alpha*del >= 0
        n = size(x,1);
       % n_0 = size(del,1)
        alpha = b-a/2;
        while 1
            h = 0;
            for i = 1:n
                if alpha*del(i) < -x(i)
                    h = 1;
                    break
                end
            end
            if h == 1
                b = alpha;
            else
                a = alpha;
            end
            alpha = (b-a)/2;
            if b-a < tol
                break;
            end
        end    
end
