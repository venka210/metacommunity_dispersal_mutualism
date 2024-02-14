clear all

syms x r_x K_x a_xy del_x y r_y K_y a_yx del_y m a q del_m d_m


eqns = [r_x*((K_x-x-a_xy*y)/K_x) - a*m - del_x == 0, r_y*((K_y-y-a_yx*x)/K_y) - del_y == 0, a*q*x - del_m - d_m*m == 0];
S = solve(eqns,[x, y,m]);
num_sol = solve(subs(eqns,[r_x K_x a_xy del_x r_y K_y a_yx del_y a d_m], [5,200,0.73,0.01,5,200,0.60,0.03,1,1]),[x,y,m]); 

