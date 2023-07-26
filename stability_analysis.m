clear all

%% parameter definition
syms p_x p_y p_m c_x c_y c_m e_x e_y e_m lambda mu

%% equilibrium expressions

eqns = [((c_x*p_x+lambda*p_m)*(1-p_x))-(e_x*p_x)-(mu*p_m) == 0, c_y-c_y*p_y-e_y == 0, (c_m*p_m*(p_x-p_m))-(e_m*p_m) == 0];
S = solve(eqns,[p_x,p_y,p_m])
