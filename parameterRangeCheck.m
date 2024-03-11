%clear all
function [valid_parameters, invalid_indices, good_bad_vals] = parameterRangeCheck(deltaValsRange,qValsRange,fVal)
% Set up a grid of parameter values
%delta_values = linspace(deltaValsRange(1), deltaValsRange(2), 24); % Adjust the number of points as needed
delta_values = deltaValsRange;
q_values = qValsRange;
%q_values = linspace(qValsRange(1), qValsRange(2), 20); % Adjust the number of points as needed

aux_var_vals = {5,200,0.73,0.01,5,200,0.60,0.03,1,1,fVal};
[r_x, K_x, alpha_xy, delta_x, r_y, K_y, alpha_yx, delta_y, a, d_m, f] = aux_var_vals{:};
%%
% Initialize a matrix to store valid parameter combinations
valid_parameters = [];
invalid_indices = [];
good_bad_vals = [];

% Check each combination of parameters
for i = 1:size(delta_values,2)
    for j = 1:size(q_values,2)
        current_parameters = [delta_values(i); q_values(j)];
        % Check if the parameters satisfy constraints
        if all(constraintFunction(current_parameters) <= 0) && all(eqm_density_fn(current_parameters, r_x, K_x, alpha_xy, delta_x, r_y, K_y, alpha_yx, delta_y, a, d_m, f) >= 0)
            good_bad_vals = [good_bad_vals, 1];
            valid_parameters = [valid_parameters, current_parameters]; % we might ultimately need a function/command that tells us which of the species has out of bound densities 
        else
            invalid_indices = [invalid_indices, [i;j]];
            good_bad_vals = [good_bad_vals, 0];
        end
    end
end
end

% disp('Valid parameter combinations:');
% disp(valid_parameters);
