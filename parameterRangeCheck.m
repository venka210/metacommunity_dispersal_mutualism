clear all
% Set up a grid of parameter values
delta_values = linspace(0, 15, 100); % Adjust the number of points as needed
q_values = linspace(0.0, 1.0, 100); % Adjust the number of points as needed

aux_var_vals = {5,200,0.73,0.01,5,200,0.60,0.03,1,1,0.7};
[r_x K_x alpha_xy delta_x r_y K_y alpha_yx delta_y a d_m f] = aux_var_vals{:};
%%
% Initialize a matrix to store valid parameter combinations
valid_parameters = [];

% Check each combination of parameters
for param1 = delta_values
    for param2 = q_values
        current_parameters = [param1; param2];
        
        % Check if the parameters satisfy constraints
        if all(constraintFunction(current_parameters) <= 0) && all(eqm_density_fn(current_parameters, r_x, K_x, alpha_xy, delta_x, r_y, K_y, alpha_yx, delta_y, a, d_m, f) >= 0)
            valid_parameters = [valid_parameters, current_parameters];
        end
    end
end

disp('Valid parameter combinations:');
disp(valid_parameters);