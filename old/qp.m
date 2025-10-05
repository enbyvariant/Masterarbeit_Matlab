function [A, b, c, H, xl, xu] = qp()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
A = [1 -2 -1 0 0; -1 -2 0 -1 0; -1 2 0 0 -1];
b = [-2; -6; -2];
c = [-2; -5; 0; 0; 0];
H = [2 0 0 0 0; 0 2 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0];
xl = zeros(5,1);
xu = 10000*ones(5,1);
end
