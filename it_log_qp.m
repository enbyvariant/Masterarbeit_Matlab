function [data] = it_log_qp(iter, it, nlp, sigma, alpha, help, dim)
% Iteration Log 
en = ones(dim.n,1);
er = ones(dim.r,1);

data.iter = iter;
%data.dual_gap = nlp.c'*it.x - (nlp.b'*it.y +nlp.xl'*it.wl -nlp.xu'*it.wu + nlp.cl'*it.zl -nlp.cu'*it.zu);
pri_fea = [nlp.A*it.x-nlp.b; it.bound_cl*(nlp.C*it.x - nlp.cl - it.tl); ...
    it.bound_cu*(-nlp.C*it.x + nlp.cu - it.tu); it.bound_xl*(it.x - nlp.xl - it.sl);...
    it.bound_xu*(-it.x + nlp.xu - it.su)];
data.pri_fea_A = full(norm(nlp.A*it.x-nlp.b));
data.pri_fea = full(norm(pri_fea));
data.dual_fea = full(norm(nlp.H*it.x + nlp.c - nlp.A'*it.y - nlp.C'*(it.zl - it.zu) - (it.wl -it.wu)));
compl = [it.bound_xl*(help.Wl *it.sl - sigma*it.mu*en); it.bound_xu*(help.Wu *it.su - sigma*it.mu*en);...
    it.bound_cl*(help.Zl*it.tl - sigma*it.mu*er); it.bound_cu*(help.Zu*it.tu - sigma*it.mu*er)];
data.compl = full(norm(compl));
data.mu = it.mu;
data.alpha = full(alpha);
disp(data);

end