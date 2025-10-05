function [nlp] = iterate(it,p, nlp, dim)
   
    nlp.C = zeros(dim.r,dim.n);
    nlp.c_i = zeros(dim.r,1);
    for i = 1:dim.r
        nlp.c_i(i) = p.c_i(i).eval_at(it.x);
        nlp.C(i,:) = p.c_i(i).gradient_at(it.x)';
    end
    nlp.A = zeros(dim.m,dim.n);
    nlp.c_e = zeros(dim.m,1);
    for i = 1:dim.m
        nlp.c_e(i) = p.c_e.eval_at(it.x);
        nlp.A(i,:) = p.c_e(i).gradient_at(it.x)';
    end
    nlp.H = p.f.hessian_at(it.x);
    nlp.obj = p.f.eval_at(it.x);
    nlp.grad = p.f.gradient_at(it.x);


end
