function [A, b, H, c, c_0, xl, xu] = get_input()
    
    prob = cutest_setup();
    n = prob.n;
    m = prob.m;
    b = sparse(-cutest_cons(zeros(n,1)));
    xl = sparse(prob.bl);
    xu = sparse(prob.bu);
    H = sparse(cutest_hess(prob.x, prob.v));
    c_0 = cutest_obj(zeros(n,1));
    A = zeros(m,n);
    c = zeros(n,1);
    for i = 1:n
        A(1:m,i) = cutest_cons([zeros(i-1,1); 1; zeros(n-i,1)])+b;
        c(i) = cutest_obj([zeros(i-1,1); 1; zeros(n-i,1)])- 1/2*H(i,i)-c_0;
    end
    A = sparse(A);
    c = sparse(c);
    % x = zeros(2*n,1);
    % for i = 1:n
    %     if prob.x(i) < 0
    %         x(i) = 0;
    %         x(n+i) = -prob.x(i);
    %     else
    %         x(i) = prob.x(i);
    %         x(n+i) = 0;
    %     end
    % end
    % 
    % if all(prob.bu > 10^5) && all(prob.bl < -10^5)
    %     A = [A -A];
    %     b = -b;
    %     c = [c; -c];
    %     H = [H -H; -H H];
    %     lambda = prob.v;
    %     s = rand(2*n,1);
    % elseif all(prob.bl == 0) && all(prob.bu > 10^5)
    %     x = rand(n,1);
    %     b = -b;
    %     lambda = rand(m,1);
    %     s = rand(n,1);
    % elseif all(prob.bl >= 0)
    %     x = rand(n,1);
    %     b = -b;
    %     lambda = rand(m,1);
    %     s = rand(n,1);
    % 
    %     % A = [A zeros(m,2*n); eye(n) eye(n) zeros(n); eye(n) zeros(n) -eye(n)];
    %     % b = [-b; prob.bu; prob.bl];
    %     % c = [c; zeros(2*n,1)];
    %     % H = [H zeros(n,2*n); zeros(2*n,3*n)];
    %     % x = rand(3*n,1);
    %     % lambda = rand(m+ 2*n,1);
    %     % s = rand(3*n,1);
    % else
    %     A = [A -A zeros(m,n) zeros(m,n); eye(n) -eye(n) eye(n) zeros(n); -eye(n) eye(n) zeros(n) eye(n)];
    %     b = [-b; prob.bu; -prob.bl];
    %     c = [c; -c; zeros(2*n,1)];
    % 
    %     H = [H -H zeros(n, 2*n); -H H zeros(n,2*n); zeros(2*n,4*n)];
    %     y = prob.bu - prob.x;
    %     z = prob.x - prob.bl;
    %     x = [x; y; z];
    %     lambda = [prob.v; rand(2*n,1)];
    %     s = rand(4*n,1);
    % end

    cutest_terminate;
end