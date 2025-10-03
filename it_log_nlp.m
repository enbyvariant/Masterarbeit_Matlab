function [data] = it_log_nlp(iter, it, nlp, sigma, alpha, help, dim,p)
% Iteration Log 
en = ones(dim.n,1);
er = ones(dim.r,1);

data.iter = iter;
pri_fea = [nlp.c_e; help.rho_l; help.rho_u; ...
        help.beta_l; help.beta_u];
data.pri_fea = full(norm(pri_fea));
dual_fea = nlp.grad - nlp.A'*it.y - nlp.C'*(it.zl-it.zu) - (it.bound_xl*it.wl - it.bound_xu*it.wu);
data.dual_fea = full(norm(dual_fea));
compl = [it.bound_xl*(help.Wl *it.sl); it.bound_xu*(help.Wu *it.su);...
    it.bound_cl*(help.Zl*it.tl); it.bound_cu*(help.Zu*it.tu)];
compl_p = [it.bound_xl*(help.Wl *it.sl - sigma*it.mu*en); it.bound_xu*(help.Wu *it.su - sigma*it.mu*en);...
    it.bound_cl*(help.Zl*it.tl - sigma*it.mu*er); it.bound_cu*(help.Zu*it.tu - sigma*it.mu*er)];
data.compl = full(norm(compl));
data.compl_p = full(norm(compl_p));
data.mu = it.mu;
data.alpha = full(alpha);
data.obj = nlp.obj;
disp(data);

end