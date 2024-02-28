function [c, ceq] = constraintFunction(parameters)
    % Lower and upper bounds of parameters
    lb = [0; 0];           % Lower bounds
    ub = [100; 1];         % Upper bounds

    % Inequality constraints (c <= 0)
    c = [-parameters(1) + lb(1);     % Parameter 1 must be >= 0
         -ub(1) + parameters(1);     % Parameter 1 must be <= 1
         -parameters(2) + lb(2);     % Parameter 2 must be >= 0
         -ub(2) + parameters(2)];    % Parameter 2 must be <= 100

    % No equality constraints
    ceq = [];
end
