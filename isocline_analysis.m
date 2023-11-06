%%%%%%%%%%%% Isocline analysis of the local patch dynamics model for
%%%%%%%%%%%% generalist or specialist frugivore
clearvars
%% Define the parameters
syms x y m r_x r_y alpha_xy alpha_yx K_x K_y del_x del_y d_m a del_m q f

%% Isocline equations
% Define the differential equations
dx = r_x.*x*((K_x-x-alpha_xy.*y)./K_x)-del_x.*x-a.*m*f.*x;
dy = r_y.*y.*((K_y-y-alpha_yx.*x)./K_y)- del_y.*y - a.*m.*(1-f).*y;
dm = a*q.*(f.*x+(1-f).*y) - del_m.*m-d_m.*m;

% Solve for the isocline equations
isocline_x = simplify(solve(dx==0, x));
isocline_y = simplify(solve(dy==0, y));
isocline_m = simplify(solve(dm==0, m));


%%
r_x = 5; r_y = 5;

alpha_xy = 0.73; alpha_yx = 0.60; %really only affects y's local density and dispersal

K_x = 200; K_y = 200; %no point touching this

del_x = 0.01; del_y = 0.05;
d_m = 0.3; 
a = 2.0; %capture rate of plant by frugivore
del_m = 10; %I'm only not starting from zero because the computational costs are absurd
q = 0.3;
f = 0.6;


%dxnewreg = matlabFunction(dx); dynewreg = matlabFunction(dy); dmnewreg = matlabFunction(dm);

%% Define the grid
[x, y, m] = meshgrid(0:20:K_x,0:20:K_y,0:2:K_x/10);
dxnew = double(simplify(subs(dx)));dynew = double(simplify(subs(dy)));dmnew = simplify(subs(dm));

%% Create isoclines
figure;
quiver3(x, y, m, dxnew, dynew, dmnew);
axis equal
%contour(x, y,dmnew, 'r');
xlabel('x');
ylabel('y');
title('Isoclines for Coupled Differential Equations');
legend('Direction Field', 'Location', 'NorthWest');
hold off;