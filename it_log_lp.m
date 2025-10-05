function [data] = it_log_lp(iter, it, nlp, sigma, alpha_pri, alpha_dual , help, dim)
% Iteration Log for LPs
% Returns objective, primal and dual step lengths, KKT error terms, 
% duality measure, duality gap and current x value for dim 2

en = ones(dim.n,1);
er = ones(dim.r,1);

% current iteration number
data.iter = iter;
% calculate duality gap
data.dual_gap = nlp.c'*it.x - (nlp.b'*it.y +nlp.xl'*it.wl -nlp.xu'*it.wu + nlp.cl'*it.zl -nlp.cu'*it.zu);
% current constraint error
pri_fea = [nlp.A*it.x-nlp.b; help.rho_l; ...
    help.rho_u; help.beta_l; help.beta_u];
data.pri_fea = norm(pri_fea);
% current stationarity error
data.dual_fea = norm(nlp.c - nlp.A'*it.y - nlp.C'*(it.zl - it.zu) - (it.wl -it.wu));
% current complementarity error (perturbed and unperturbed)
compl_p = [it.bound_xl*(help.Wl *it.sl - sigma*it.mu*en); it.bound_xu*(help.Wu *it.su - sigma*it.mu*en);...
    it.bound_cl*(help.Zl*it.tl - sigma*it.mu*er); it.bound_cu*(help.Zu*it.tu - sigma*it.mu*er)];
compl = [it.bound_xl*(help.Wl *it.sl); it.bound_xu*(help.Wu *it.su);...
    it.bound_cl*(help.Zl*it.tl); it.bound_cu*(help.Zu*it.tu)];
data.compl_p = full(norm(compl_p));
data.compl = full(norm(compl));
% current duality measure
data.mu = it.mu;
% current step lengths
data.alpha_p = alpha_pri;
data.alpha_d = alpha_dual;
% current iterate for dim=2
if dim.n == 2
    data.x1 = it.x(1);
    data.x2 = it.x(2);
end
disp(data);
end