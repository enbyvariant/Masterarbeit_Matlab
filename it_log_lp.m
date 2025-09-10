function [data] = it_log_lp(iter, it, nlp, sigma, alpha_pri, alpha_dual , help, dim)
% Iteration Log 
en = ones(dim.n,1);
er = ones(dim.r,1);

data.iter = iter;
data.dual_gap = nlp.c'*it.x - (nlp.b'*it.y +nlp.xl'*it.wl -nlp.xu'*it.wu + nlp.cl'*it.zl -nlp.cu'*it.zu);
pri_fea = [nlp.A*it.x-nlp.b; help.rho_l; ...
    help.rho_u; help.beta_l; help.beta_u];
data.pri_fea = norm(pri_fea);
data.dual_fea = norm(nlp.c - nlp.A'*it.y - nlp.C'*(it.zl - it.zu) - (it.wl -it.wu));
compl = [help.Wl *it.sl - sigma*it.mu*en; help.Wu *it.su - sigma*it.mu*en;...
    help.Zl*it.tl - sigma*it.mu*er; help.Zu*it.tu - sigma*it.mu*er];
data.mu = it.mu;
data.compl = norm(compl);
data.alpha_p = alpha_pri;
data.alpha_d = alpha_dual;
disp(data);

end