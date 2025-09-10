function [p] = test_problem()
% create test problem with symbolic
%   variables.

p.n = 2;
p.m = 2;

p.X = sym("x", [2,1]);
f(p.X) = sin(p.X(1)) + p.X(2)^2;
p.f = Symbolic_Function(f,p.X);

p.c_e = [];
c_i = [-p.X(1)-p.X(2) + 2;-p.X(1) + p.X(2) + 1];
p.c_i = [Symbolic_Function(c_i(1),p.X);
            Symbolic_Function(c_i(2),p.X)];
p.bl = zeros(2,1);
p.bu = 10000*ones(2,1);

p.cl = [0;0];
p.cu = [10000;10000];
end