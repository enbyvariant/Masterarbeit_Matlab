function [gamma] = matrix_factors(help, gamma, dim, nlp)
% Compute factors for 
%   Regularization and Convexification of Omega

if isempty(gamma)
    gamma.con = 0;
    gamma.reg = 0;
    gamma.b = 0.8;
    gamma.e = 0.1;
end
Hphi = nlp.H + help.PHI+ gamma.reg*eye(dim.n);

mat = [Hphi   nlp.A'  nlp.C';
       nlp.A   zeros(dim.m)    zeros(dim.m,dim.r);
       nlp.C  zeros(dim.r,dim.m)   -help.PSI1];

d = eig(mat);
n_plus = 0;
n_minus = 0;
n_0 = 0;
for i = 1:length(d)
    if d(i) > 0
        n_plus = n_plus + 1;
    elseif d(i) < 0
        n_minus = n_minus + 1;
    else
        n_0 = n_0 + 1;
    end
end
d_minus = zeros(n_minus,1);
i = 1;
for e = d
    if e < 0
        d_minus(i) = e;
        i = i + 1;
    end
end
        
d_minus = sort(d_minus);
if n_plus < dim.n
    disp('blub')
    if dim.n - n_plus == n_0
        gamma.reg = 0.05;
    else
        gamma.reg = -1.1*d_minus(n_minus+1 -(dim.n - n_plus - n_0));
    end
else
    gamma.reg = 0;
end
          

% [L,~] = lu(nlp.A);
% flag = find(all(L==0));
% r = rank(full(nlp.A));
% if r < dim.m
%     disp("con")
%     gamma.con = 10^(-8) * gamma.e * it.mu^(gamma.b);
% else
%     gamma.con = 0;
% end


end