classdef Symbolic_Function
    %SYMBOLIC_FUNCTION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        f
        grad_f
        hessian_f
        symv
        dim
    end
    
    methods
        function obj = Symbolic_Function(f,symv)
            obj.f = matlabFunction(f, 'Vars', {symv});
            obj.symv = symv;
            obj.dim = length(argnames(f));
            %obj.grad_f = gradient(obj.f,obj.symv);
            obj.grad_f = gradient(f,obj.symv);
            obj.grad_f = matlabFunction(obj.grad_f, 'Vars', {symv});
            obj.hessian_f = hessian(f,obj.symv);
            obj.hessian_f = matlabFunction(obj.hessian_f, 'Vars', {symv});

        end
        
        function column_vector = assert_valid_input(~,x)
            s = size(x);
            assert(s(1) == 1 || s(2) == 1, "Input cannot be a matrix!");
            if (s(1) == 1) 
                column_vector = x;
            else
                column_vector = x';
            end
        end

        function val = eval_at(obj,x)
            % Evaluates f (R^n -> R) at a point x
            % Takes x as column or row vector
            %x = obj.assert_valid_input(x);
            %val = double(subs(obj.f,obj.symv,x));
            val = obj.f(x);
        end

        function g = gradient_at(obj,x)
        % Evaluates gradient at a point x
        % Takes x as column or row vector
            %x = obj.assert_valid_input(x);
            %g = double(subs(obj.grad_f,obj.symv,x));
            g = obj.grad_f(x);
        end
        
        function H = hessian_at(obj,x)
            %x = obj.assert_valid_input(x);
            %H = double(subs(hessian(obj.f,obj.symv),obj.symv,x));
            H = obj.hessian_f(x);
        end
    end
end

