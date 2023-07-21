function dydt = vectorized_BetweenPatchDynamics(t,y, c_x, c_y, c_m, e_x, e_y, e_m, mu, lambda,n)
% y(1,:) = occupied patch fraction by species x (p_x)
% y(2,:) = occupied patch fraction by species y (p_y)
% y(3,:) = occupied patch fraction by species m (p_y)
dydt = zeros(3,n);
y = reshape(y,[],n);
dydt(1,:) = (c_x.*y(1,:)+lambda.*y(3,:)).*(1-y(1,:))-e_x.*y(1,:)-mu.*y(3,:);%removing the c_m component
dydt(2,:) = c_y.*y(2,:).*(1-y(2,:))-e_y.*y(2,:);
dydt(3,:) = c_m.*y(3,:).*(y(1,:)-y(3,:))-e_m.*y(3,:);%-(e_x+mu).*y(3,:);
dydt = dydt(:);
end
