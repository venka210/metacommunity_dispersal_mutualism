function stateVariables = eqm_density_fn(parameters, r_x, K_x, alpha_xy, delta_x, r_y, K_y, alpha_yx, delta_y, a, d_m, f)
stateVariable1 = (K_x*d_m*r_x*r_y - K_x*d_m*delta_x*r_y - K_y*a*alpha_xy*parameters(1)*r_x + K_y*alpha_xy*d_m*delta_y*r_x + K_x*a*parameters(1)*f*r_y - K_y*alpha_xy*d_m*r_x*r_y - K_x*K_y*a^2*delta_x*parameters(2) + K_x*K_y*a^2*parameters(2)*r_x - K_x*K_y*a^2*delta_x*f^2*parameters(2) - K_x*K_y*a^2*delta_y*f^2*parameters(2) + K_x*K_y*a^2*f^2*parameters(2)*r_x + K_x*K_y*a^2*f^2*parameters(2)*r_y + K_y*a*alpha_xy*parameters(1)*f*r_x + 2*K_x*K_y*a^2*delta_x*f*parameters(2) + K_x*K_y*a^2*delta_y*f*parameters(2) - 2*K_x*K_y*a^2*f*parameters(2)*r_x - K_x*K_y*a^2*f*parameters(2)*r_y)/(d_m*r_x*r_y + K_y*a^2*parameters(2)*r_x + K_x*a^2*f^2*parameters(2)*r_y + K_y*a^2*f^2*parameters(2)*r_x - alpha_xy*alpha_yx*d_m*r_x*r_y - 2*K_y*a^2*f*parameters(2)*r_x - K_x*a^2*alpha_yx*f*parameters(2)*r_y - K_y*a^2*alpha_xy*f*parameters(2)*r_x + K_x*a^2*alpha_yx*f^2*parameters(2)*r_y + K_y*a^2*alpha_xy*f^2*parameters(2)*r_x);
stateVariable2 = -(K_y*d_m*delta_y*r_x - K_y*a*parameters(1)*r_x - K_y*d_m*r_x*r_y - K_x*alpha_yx*d_m*delta_x*r_y + K_y*a*parameters(1)*f*r_x + K_x*alpha_yx*d_m*r_x*r_y + K_x*K_y*a^2*delta_x*f^2*parameters(2) + K_x*K_y*a^2*delta_y*f^2*parameters(2) - K_x*K_y*a^2*f^2*parameters(2)*r_x - K_x*K_y*a^2*f^2*parameters(2)*r_y + K_x*a*alpha_yx*parameters(1)*f*r_y - K_x*K_y*a^2*delta_x*f*parameters(2) + K_x*K_y*a^2*f*parameters(2)*r_x)/(d_m*r_x*r_y + K_y*a^2*parameters(2)*r_x + K_x*a^2*f^2*parameters(2)*r_y + K_y*a^2*f^2*parameters(2)*r_x - alpha_xy*alpha_yx*d_m*r_x*r_y - 2*K_y*a^2*f*parameters(2)*r_x - K_x*a^2*alpha_yx*f*parameters(2)*r_y - K_y*a^2*alpha_xy*f*parameters(2)*r_x + K_x*a^2*alpha_yx*f^2*parameters(2)*r_y + K_y*a^2*alpha_xy*f^2*parameters(2)*r_x);
stateVariable3 = -(parameters(1)*r_x*r_y + K_y*a*delta_y*parameters(2)*r_x - K_y*a*parameters(2)*r_x*r_y - alpha_xy*alpha_yx*parameters(1)*r_x*r_y - K_x*a*alpha_yx*delta_x*parameters(2)*r_y + K_x*a*delta_x*f*parameters(2)*r_y - K_y*a*delta_y*f*parameters(2)*r_x + K_x*a*alpha_yx*parameters(2)*r_x*r_y - K_x*a*f*parameters(2)*r_x*r_y + K_y*a*f*parameters(2)*r_x*r_y + K_x*a*alpha_yx*delta_x*f*parameters(2)*r_y - K_y*a*alpha_xy*delta_y*f*parameters(2)*r_x - K_x*a*alpha_yx*f*parameters(2)*r_x*r_y + K_y*a*alpha_xy*f*parameters(2)*r_x*r_y)/(d_m*r_x*r_y + K_y*a^2*parameters(2)*r_x + K_x*a^2*f^2*parameters(2)*r_y + K_y*a^2*f^2*parameters(2)*r_x - alpha_xy*alpha_yx*d_m*r_x*r_y - 2*K_y*a^2*f*parameters(2)*r_x - K_x*a^2*alpha_yx*f*parameters(2)*r_y - K_y*a^2*alpha_xy*f*parameters(2)*r_x + K_x*a^2*alpha_yx*f^2*parameters(2)*r_y + K_y*a^2*alpha_xy*f^2*parameters(2)*r_x);
stateVariables = [stateVariable1; stateVariable2; stateVariable3]
end
