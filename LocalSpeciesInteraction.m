function dydt = LocalSpeciesInteraction(t,y, r_x,r_y, alpha_xy, alpha_yx, K_x, K_y, del_x, del_y, del_m, a, q, d_m)
% y(1) = population of species x
% y(2) = population of species y
% y(3) = population of species z

% system of ODEs
dydt = zeros(3,1);
dydt(1) = r_x.*y(1).*((K_x - y(1) - alpha_xy.*y(2))/(K_x)) - a.*y(1).*y(3) - del_x.*y(1);
dydt(2) = r_y.*y(2).*((K_y - y(2) - alpha_yx.*y(1))/(K_y)) - del_y.*y(2);
dydt(3) = a.*q.*y(3).*y(1) - del_m.*y(3) - d_m.*y(3);
end