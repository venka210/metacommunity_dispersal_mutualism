%%%%%%%% This code vectorizes the ODEs that need to be solved for multiple
%%%%%%%% varying parameters (here 'a' and 'del_m')
clear vars

%% parameter definition
r_x = 5; r_y = 5;

alpha_xy = 0.73; alpha_yx = 0.60; %really only affects y's local density and dispersal

K_x = 200; K_y = 200; %no point touching this

del_x = 0.01; del_y = 0.05;
q = 0.55; d_m = 0.3; 

del_m = 3:0.2:22; %I'm only not starting from zero because the computational costs are absurd
a = 0.46:0.01:0.81;%0.7-1.1 seems to work for this fig.

[Del,A] = meshgrid(del_m,a);

Del_col = Del(:); A_col = A(:);

num_combinations = numel(A_col); 

k_eff = 0.1; %efficiency of dispersing seeds to habitable patches

z_x = 0.7; z_y = 0.7; z_m = 0.4; %scaling factors for patch extinction rates. changing z_m relative to z_x and z_m does not change qual. change results

e_xmin = 0.05;e_ymin = 0.05; e_mmin = 0.05; e_mxmin = e_xmin;

tspan = [0,1000];

k_x = 0.08; k_y = 0.08; k_m = 0.08;

x_init = 0.1; y_init = 0.1; m_init = 0.1;
    
spp_init_no_m = [x_init y_init 0];
spp_init = [x_init y_init m_init];


%variable collectors across parameter sweeps

occupancy_del_m = zeros(size(del_m,1),3);
eta_check = zeros(size(del_m,1),1);
mu_collector = zeros(size(del_m,1),1);
gamma = zeros(size(del_m,1),1);
cm_collector = zeros(size(del_m,1),1);
eta2 = zeros(size(del_m,1),1);
em_collector = zeros(size(del_m,1),1);
lambda_collector = zeros(size(del_m,1),1);
%%

options = odeset('NonNegative',[1,2,3]);
[t_patch_no_m,local_dens_no_m] = ode45(@(t,y) vectorized_LocalSpeciesInteraction(t,y,r_x,r_y,alpha_xy,alpha_yx, K_x, K_y, del_x, del_y, Del_col', A_col', q, d_m, num_combinations), tspan/10, repmat(spp_init,1,num_combinations), options);
[t_patch,local_dens] = ode45(@(t,y)vectorized_LocalSpeciesInteraction(t,y,r_x,r_y,alpha_xy,alpha_yx, K_x, K_y, del_x, del_y,Del_col', A_col', q, d_m, num_combinations), tspan, repmat(spp_init, 1, num_combinations), options);

num_timepoints_no_m = length(t_patch_no_m);
num_timepoints = length(t_patch);
 
local_dens_no_m_3d = reshape(local_dens_no_m(end,:),[],num_combinations);
local_dens_3d = reshape(local_dens(end,:),[],num_combinations);




%      figure()
%      plot(t_patch,local_dens)
%% metacommunity dynamics
e_x = e_xmin.*((K_x./(local_dens_no_m_3d(1,:))).^z_x);e_y = e_ymin.*(K_y./(local_dens_3d(2,:))).^z_y;
e_mx = e_xmin.*(K_x./(local_dens_3d(1,:))).^z_x;
e_m = e_mmin.*(K_x./(local_dens_3d(3,:))).^z_m; %assume max population size of mutualist is similar to that of x and y


mu = e_mx - e_x;

tspan_meta = [0, 10000];

c_x = (k_x*del_x).*local_dens_no_m_3d(1,:);
c_y = k_y*del_y.*local_dens_3d(2,:);
c_m = k_m*Del_col'.*local_dens_3d(3,:);%repelem(del_m,numel(a)) can also be used in place of Del_col' if u feel funky
c_mx =(k_x*del_x+k_m*k_eff.*A_col'*(1-q).*Del_col'.*local_dens_3d(3,:)).*local_dens_3d(1,:);

lambda = c_mx-c_x;

options = odeset('NonNegative',1);
[t_syst, frac_occup] = ode45(@(t,y)vectorized_BetweenPatchDynamics(t,y, c_x, c_y, c_m, e_x, e_y, e_m, mu, lambda, num_combinations), tspan_meta, local_dens_3d(:),options);


num_tpts_t_syst = length(t_syst); frac_occup_3d = reshape(frac_occup(end,:), [],num_combinations);
% 
% eta_check = 0.5.*sqrt((((lambda.*(1+(e_m./c_m))+c_x)-(e_x+mu))./(c_x+lambda)).^2 + 4.*(e_m./c_m).*((mu-lambda)./(c_x+lambda)));
%      figure()
%      plot(t_syst(end-100:end),frac_occup(end-100:end,1:2));

% mu_collector(i,1) = mu;
% cm_collector(i,1) = c_m;
% em_collector(i,1) = e_m;
% lambda_collector(i,1) = lambda;
% occupancy_del_m(i, :) = frac_occup(end,:);
% gamma(i,1) = (0.5.*((lambda*(1+(e_m./c_m))+c_x)-(e_x+mu)))./(c_x+lambda);%0.5*(1-((mu-e_x-e_m)/(c_x+c_m)));
% eta2(i,:) = 4*(e_m./c_m)*((mu-lambda)./(c_x+lambda));
save ("vectorized_a_del_m_varied.mat")
%% surface plot
figure()
surf(Del, A, reshape(frac_occup_3d(1,:),numel(a),numel(del_m)))
colormap parula; freezeColors;
hold on
surf(Del, A, reshape(frac_occup_3d(2,:),numel(a),numel(del_m)))
colormap winter; freezeColors;
%surf(Del, A,occupancy_del_m(:,:,1))
xlabel('mutualist dispersal rate')
ylabel ('predation rate')
%xlim([1.0 29.0]);
zlabel('fraction of patches occupied')
%xline([1.0 2.6, 10.4, 26.9],'--',{'Exploitative (x extinct)','Mutualism (y fitter)','Mutualism (x fitter)', 'Mutualism (y fitter)'})
title('Fraction of patches occupied vs mutualist dispersal rate and predation rate')
legend('Species with mutualist (x)')%, 'Species without mutualist (y)', 'mutualist (m)', 'location', 'best' )
%fig1name = sprintf('occupancy_vs_del_m.jpeg');
print('vectorized_occupancy_3d_vs_del_m','-djpeg','-r600')
