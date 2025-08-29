function [gamma, mat] = matrix_factors(help, gamma, dim, nlp,it)
% Compute factors for 
%   Regularization and Convexification of Omega

mat = nlp.H + help.PHI;
[~,flag] = chol(mat);
if flag
    gamma.reg = gamma.reg/2;
    if gamma.reg == 0
        gamma.reg = 10^(-4);
    end
    while true
        mat = nlp.H + help.PHI + gamma.reg*eye(dim.n);
        [~,flag] = chol(mat);
        if not(flag)
            break;
        end
        gamma.reg = 10*gamma.reg;
    end
else
    gamma.reg = 0;
end

% [L,~] = lu(nlp.A);
% flag = find(all(L==0));
if rank(full(nlp.A)) < dim.m
    gamma.con = 10^(-8) * gamma.e * it.mu^(gamma.b);
else
    gamma.con = 0;
end


end