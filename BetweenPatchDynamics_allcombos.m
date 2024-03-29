function dydt = BetweenPatchDynamics_allcombos(t,y, c_x0, c_xm, c_xym, c_xy, c_y0, c_ym, c_yxm, c_yx, c_mx, c_mxy, c_my, e_x0, e_xy, e_xm, e_xym, e_y0, e_yx, e_ym, e_yxm, ...
    e_mx, e_mxy, e_my)
% y(1) = occupied patch fraction by species x
% y(2) = occupied patch fraction by species y
% y(3) = occupied patch fraction by species m
dydt = zeros(3,1);

dydt(1) = y(1)*((c_x0*(1-y(1))-e_x0)*(1-(y(3)+y(2))+y(2)*y(3)) + (c_xm*(1-y(1))-e_xm)*y(3)*(1-y(2)) + (c_xy*(1-y(1))-e_xy)*y(2)*(1-y(3)) + (c_xym*(1-y(1))-e_xym)*y(2)*y(3));%removing the c_m component
dydt(2) = y(2)*((c_y0*(1-y(2))-e_y0)*(1-(y(3)+y(1))+y(1)*y(3)) + (c_ym*(1-y(2))-e_ym)*y(3)*(1-y(1)) + (c_yx*(1-y(2))-e_yx)*y(1)*(1-y(3)) + (c_yxm*(1-y(2))-e_yxm)*y(1)*y(3));
dydt(3) = y(3)*((c_mx*(1-y(3))-e_mx)*(1-y(2))*y(1)+ (c_my*(1-y(3))-e_my)*(1-y(1))*y(2)+ (c_mxy*(1-y(3))-e_mxy)*y(2)*y(1));%-(e_x+mu)*y(3);

end

% y(1,:) = occupied patch fraction by species x
% y(2,:) = occupied patch fraction by species y
% y(3,:) = occupied patch fraction by species m
% dydt = zeros(3,n);
% y = reshape(y,[],n);
% dydt(1,:) = y(1,:).*((c_x0.*(1-y(1,:))-e_x0).*(1-(y(3,:)+y(2,:))+y(2,:).*y(3,:)) + (c_xm.*(1-y(1,:))-e_xm).*y(3,:).*(1-y(2,:))) + (c_xy.*(1-y(1,:))-e_xy).*y(2,:).*(1-y(3,:)) + (c_xym.*(1-y(1,:))-e_xym).*y(2,:).*y(3,:);%removing the c_m component
% dydt(2,:) = y(2,:).*((c_y0.*(1-y(2,:))-e_y0).*(1-(y(3,:)+y(1,:))+y(1,:).*y(3,:)) + (c_ym.*(1-y(2,:))-e_ym).*(y(3,:).*(1-y(1,:))) + (c_yx.*(1-y(2,:))-e_yx).*(y(1,:).*(1-y(3,:))) + (c_yxm.*(1-y(2,:))-e_yxm).*(y(1,:).*y(3,:)));
% %dydt(3,:) = y(3,:).*((1-y(2,:)).*(c_mx.*(y(1,:)-y(3,:))-e_mx)+y(2,:).*(c_mxy.*(y(1,:)-y(3,:))-e_mxy));
% dydt(3,:) = y(3,:).*((c_mx.*(1-y(3,:))-e_mx).*(1-y(2,:)).*y(1,:) + (c_my.*(1-y(3,:))-e_my).*(1-y(1,:)).*y(2,:) + (c_mxy.*(1-y(3,:))-e_mxy).*y(1,:).*y(2,:));
% dydt = dydt(:);
