function dydt = BetweenPatchDynamics_pillai(t,y, c_x, c_y, c_m, e_x, e_y, e_m, mu)
% y(1) = occupied patch fraction by species x
% y(2) = occupied patch fraction by species y
% y(3) = occupied patch fraction by species m
dydt = zeros(3,1);
dydt(1) = (c_x*y(1) + c_m*y(3)).*(1-y(1))-e_x*y(1)-mu*y(3);
dydt(2) = c_y*y(2).*(1-y(2))-e_y.*y(2);
dydt(3) = c_m*y(3).*(y(1)-y(3))-e_m*y(3)-(e_x+mu)*y(3);
end
