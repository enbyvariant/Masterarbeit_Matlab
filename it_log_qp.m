function [data] = it_log_qp(iter, it, nlp, sigma, alpha, help, dim)
% Iteration Log for QPs
% Returns objective, step length, KKT error terms, duality measure, 
% duality gap and current x value for dim 2

en = ones(dim.n,1);
er = ones(dim.r,1);

% current iteration
data.iter = iter;
% current objective
pri_obj = 1/2*it.x'*nlp.H*it.x + nlp.c'*it.x;
data.obj = pri_obj;
% compute dual objective and duality gap
x_inf = lsqminnorm(nlp.H, -nlp.c + nlp.A'*it.y + nlp.C'*(it.zl - it.zu) + (it.wl - it.wu));
dual_obj = -1/2*(x_inf)'*nlp.H * (x_inf) + nlp.b'*it.y +nlp.xl'*it.wl -nlp.xu'*it.wu + nlp.cl'*it.zl -nlp.cu'*it.zu;
data.dual_gap = pri_obj - dual_obj;
% current constraint error
pri_fea = [nlp.A*it.x-nlp.b; it.bound_cl*(nlp.C*it.x - nlp.cl - it.tl); ...
    it.bound_cu*(-nlp.C*it.x + nlp.cu - it.tu); it.bound_xl*(it.x - nlp.xl - it.sl);...
    it.bound_xu*(-it.x + nlp.xu - it.su)];
data.pri_fea = full(norm(pri_fea));
% current stationarity error
data.dual_fea = full(norm(nlp.H*it.x + nlp.c - nlp.A'*it.y - nlp.C'*(it.zl - it.zu) - (it.wl -it.wu)));
% current complementarity error (perturbed and unperturbed)
compl_p = [it.bound_xl*(help.Wl *it.sl - sigma*it.mu*en); it.bound_xu*(help.Wu *it.su - sigma*it.mu*en);...
    it.bound_cl*(help.Zl*it.tl - sigma*it.mu*er); it.bound_cu*(help.Zu*it.tu - sigma*it.mu*er)];
compl = [it.bound_xl*(help.Wl *it.sl); it.bound_xu*(help.Wu *it.su);...
    it.bound_cl*(help.Zl*it.tl); it.bound_cu*(help.Zu*it.tu)];
data.compl_p = full(norm(compl_p));
data.compl = full(norm(compl));
% current duality measure
data.mu = it.mu;
% current step length
data.alpha = full(alpha);
% current iterate if problem has dim 2
if dim.n == 2
    data.x1 = it.x(1);
    data.x2 = it.x(2);
end
disp(data);

end