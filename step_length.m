function [alpha] = step_length(step, curr, eta)
%   Compute alpha such that
%   eta*curr + alpha*step = 0
%   close to maximal step length that fullfills positivity
    
    

    index = find(step < 0);
    if isempty(index)
        alpha = 1;
    else
        alpha = eta*min(curr(index)./(-step(index)));
        if alpha > 1
            alpha = 1;
        end
    end
end