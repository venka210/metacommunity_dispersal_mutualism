clear all

%% parameter definition
sym p_x p_y p_m
r_x = 5; r_y = 5;

alpha_xy = 0.73; alpha_yx = 0.60; %really only affects y's local density and dispersal

K_x = 200; K_y = 200; %no point touching this

del_x = 0.01; del_y = 0.05;
%q = 0.45;
d_m = 0.3; 
a = 2.0; %capture rate of plant by frugivore
del_m = 5; %I'm only not starting from zero because the computational costs are absurd
%a = 0.51:0.01:0.81;%0.7-1.1 seems to work for this fig.
q = 0.8;
% [Del,Q] = meshgrid(del_m,q);
% 
% Del_col = Del(:); Q_col = Q(:);

%num_combinations = numel(Q_col); 

k_eff = 1; %efficiency of dispersing seeds to habitable patches

z_x = 0.7; z_y = 0.7; z_m = 0.4; %scaling factors for patch extinction rates. changing z_m relative to z_x and z_m does not change qual. change results

e_xmin = 0.05;e_ymin = 0.05; e_mmin = 0.05; e_mxmin = e_xmin;

%tspan = [0,1000];

k_x = 0.08; k_y = 0.08; k_m = 0.08;

x_init = 0.1; y_init = 0.1; m_init = 0.1;

%% equilibrium expressions

eqn1 = 
