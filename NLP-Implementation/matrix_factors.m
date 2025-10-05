function [gamma] = matrix_factors(dim, nlp)
% Compute factors for 
%   Regularization and Convexification of Omega

% Hphi = nlp.H + help.PHI;
if dim.m ~= 0
    cell = spspaces(nlp.A', 1);
    Z = cell{1}(dim.m+1:dim.n,:);
    red_H = Z*nlp.H*Z';
else
    red_H = nlp.H;
end
[~,flag] = chol(red_H);
if flag
    d = eigs(red_H,6,"smallestreal");
    d = min(d);
    gamma = -1.1*d;
else
    gamma = 0;
end
if dim.m ~= 0
    red_H = Z*(nlp.H + gamma*eye(dim.n))*Z';
    [~,flag] = chol(red_H);
    if flag
        disp('nope');
    end
end
end
