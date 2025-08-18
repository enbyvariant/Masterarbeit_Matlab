function [it, obj, iterations, data, step] = Interior_Points_LP_alt(iter, nlp, dim)
    
    % Initial values
    n = dim.n;
    m = dim.m;
    r = dim.r;
    eta = 0.995;

    % Compute starting point
    [it, step_zero] = Init_LP(nlp, dim);
    disp(it.x)
    nlp.H = zeros(n);

    iterations = 0;
    data = zeros(iter,3);
    p = cutest_setup;


    while 1
        if iterations > iter
            break;
        end
        iterations = iterations + 1;

        % update helping variables
        [help] = helpers(dim, nlp, it);

        % solve linear equation system for the affine step
        [step_a] = LES_Solver(dim, nlp, help, it, 0, step_zero);
        %disp(step_a.del_x);
                    
        % calculate primal affine step length
        step_pri = [step_a.del_sl; step_a.del_su; step_a.del_tl; step_a.del_tu];
        curr = [it.sl; it.su; it.tl; it.tu];
        alpha_pri = step_length(step_pri,curr,1);

        % calculate dual affine step length
        step_dual = [step_a.del_wl; step_a.del_wu; step_a.del_zl; step_a.del_zu];
        curr = [it.wl; it.wu; it.zl; it.zu];
        alpha_dual = step_length(step_dual, curr, 1);
            
        % calculate affine duality measure
        sl_a = it.sl + alpha_pri*step_a.del_sl;
        su_a = it.su + alpha_pri*step_a.del_su;
        tl_a = it.tl + alpha_pri*step_a.del_tl;
        tu_a = it.tu + alpha_pri*step_a.del_tu;
        wl_a = it.wl + alpha_dual*step_a.del_wl;
        wu_a = it.wu + alpha_dual*step_a.del_wu;
        zl_a = it.zl + alpha_dual*step_a.del_zl;
        zu_a = it.zu + alpha_dual*step_a.del_zu;
        mu_aff = (sl_a'*wl_a + su_a'*wu_a + tl_a'*zl_a + tu_a'*zu_a)/(2*(n+r));
    
        % set the centering parameter sigma
        sigma = (mu_aff/it.mu)^3;
    
        [step] = LES_Solver(dim, nlp, help, it, sigma, step_a);
        %disp(step);

        % calculate primal step length
        step_pri = [step.del_sl; step.del_su; step.del_tl; step.del_tu];
        curr = [it.sl; it.su; it.tl; it.tu];
        alpha_pri = step_length(step_pri,curr,eta);

        % calculate dual step length
        step_dual = [step.del_wl; step.del_wu; step.del_zl; step.del_zu];
        curr = [it.wl; it.wu; it.zl; it.zu];
        alpha_dual = step_length(step_dual, curr, eta);

        % update iterate
        it.x = it.x + alpha_pri*step.del_x;
        it.y = it.y + alpha_dual*step.del_y;
        it.zl = it.zl + alpha_dual*step.del_zl;
        it.zu = it.zu + alpha_dual*step.del_zu;
        it.wl = it.wl + alpha_dual*step.del_wl;
        it.wu = it.wu + alpha_dual*step.del_wu;
        it.sl = it.sl + alpha_pri*step.del_sl;
        it.su = it.su + alpha_pri*step.del_su;
        it.tl = it.tl + alpha_pri*step.del_tl;
        it.tu = it.tu + alpha_pri*step.del_tu;
        disp([iterations;it.x])

        it.mu = (it.sl'*it.wl + it.su'*it.wu + it.tl'*it.zl + it.tu'*it.zu)/(2*(n+r));
        obj = cutest_obj(it.x);

        data(iterations, :) = [iterations obj it.mu];
        if it.mu < 1*10^(-15)
            break;
        end
    end
    fprintf('      Iteration  objective     mu-value  \n');
    disp(data);
    cutest_terminate
end