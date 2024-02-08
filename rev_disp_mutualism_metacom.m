clear all

%% parameter definition
r_x = 5; r_y = 5;

alpha_xy = 0.73; alpha_yx = 0.60; %really only affects y's local density and dispersal -- alpha_xy = 0.73; alpha_yx = 0.60;

K_x = 200; K_y = 200; %no point touching this

del_x = 0.01; del_y = 0.05; del_m = 0:1:15; %I'm only not starting from zero because the computational costs are absurd

a = 1; q = 1; d_m = 1; %
k_eff = 1; %efficiency of dispersing seeds to habitable patches

z_x = 0.8; z_y = 0.8; z_m = 0.5; %scaling factors for patch extinction rates. changing z_m relative to z_x and z_m does not change qual. change results

e_xmin = 0.05;e_ymin = 0.05; e_mmin = 0.05; e_mxmin = e_xmin;

tspan = [0,500];

k_x = 0.8; k_y = 0.8; k_m = 0.08;

x_init = 0.1; y_init = 0.1; m_init = 0.1;
    
spp_init_no_m = [x_init; y_init; 0];
spp_init_no_y = [x_init; 0; m_init];
spp_init_no_x = [0; y_init; m_init];
spp_init = [x_init; y_init; m_init];

%variable collectors across parameter sweeps

occupancy_del_m = zeros(size(del_m,1),3);
% eta_check = zeros(size(del_m,1),1);
% mu_collector = zeros(size(del_m,1),1);
% gamma = zeros(size(del_m,1),1);
cmx_collector = zeros(size(del_m,1),1);
% eta2 = zeros(size(del_m,1),1);
% em_collector = zeros(size(del_m,1),1);
% lambda_collector = zeros(size(del_m,1),1);


%% within patch dynamics

for i = 1:length(del_m)
    
    options = odeset('NonNegative',[1,2,3]);
    [t_patch_no_m,local_dens_no_m] = ode45(@(t,y)LocalSpeciesInteraction(t,y,r_x,r_y,alpha_xy,alpha_yx, K_x, K_y, del_x, del_y, del_m(i), a, q, d_m), tspan, spp_init_no_m, options);
    [t_patch_no_y,local_dens_no_y] = ode45(@(t,y)LocalSpeciesInteraction(t,y,r_x,r_y,alpha_xy,alpha_yx, K_x, K_y, del_x, del_y, del_m(i), a, q, d_m), tspan, spp_init_no_y,options);
    [t_patch,local_dens] = ode45(@(t,y)LocalSpeciesInteraction(t,y,r_x,r_y,alpha_xy,alpha_yx, K_x, K_y, del_x, del_y, del_m(i), a, q, d_m), tspan, spp_init, options);

%      figure()
%      plot(t_patch,local_dens)
    %% metacommunity dynamics
    tspan_meta = [0, 1000];

    e_x0 = e_xmin;e_y0 = e_ymin;
    e_xy = e_xmin*((K_x/(local_dens_no_m(end,1)))^z_x); e_yx = e_ymin*((K_y/(local_dens_no_m(end,2)))^z_y);
    e_xm = e_xmin*((K_x/(local_dens_no_y(end,1)))^z_x); e_ym = e_ymin;
    e_xym = e_xmin*(K_x/(local_dens(end,1)))^z_x; e_yxm = e_ymin*(K_y/(local_dens(end,2)))^z_y;
    e_mx = e_mmin*(K_x/(local_dens_no_y(end,3)))^z_m; e_my = 0; e_mxy = e_mmin*(K_x/(local_dens(end,3)))^z_m; %assume max population size of mutualist is similar to that of x and y
    e_my = 0;

    c_x0 = (k_x*del_x)*K_x*(1-(del_x/r_x)); c_y0 = (k_y*del_y)*K_y*(1-(del_y/r_y)); %patches with only one species
    c_xy = (k_x*del_x)*local_dens_no_m(end,1); c_yx = (k_y*del_y)*local_dens_no_m(end,2); %patches with 2 species
    c_xm = (k_x*del_x)*local_dens_no_y(end,1)+k_m*k_eff.*del_m(i).*local_dens_no_y(end,3)*(K_x - local_dens_no_y(end,1)); c_ym = c_y0; c_mx = k_m*del_m(i)*local_dens_no_y(end,3);% patches with one plant-one frugivore
    c_xym = k_x*del_x*local_dens(end,1)+k_m*k_eff*del_m(i)*local_dens(end,3)*(K_x - local_dens(end,1)); c_yxm = k_y*del_y*local_dens(end,2); c_mxy = k_m*del_m(i)*local_dens(end,3);%all species present
    c_my = 0;

    if any(e_mxy == inf) || any(e_mx ==inf)
    c_mx(find((e_mx == inf))) = 0; e_mx(find((e_mx == inf))) = 0;  
    c_mxy(find((e_mxy == inf))) = 0; e_mxy(find((e_mxy == inf))) = 0;
    end

    options2 = odeset('NonNegative',[1,2,3]);
    frac_occup_init = [0.1;0.1;0.1]; %initial patch occupancies of each species x, y, and, m respectively. 

    [t_syst, frac_occup] = ode45(@(t,y)BetweenPatchDynamics_allcombos(t,y, c_x0, c_xm, c_xym, c_xy, c_y0, c_ym, c_yxm, c_yx, c_mx, c_mxy, c_my, e_x0, e_xy, e_xm, e_xym, e_y0, e_yx, e_ym, e_yxm, ...
    e_mx, e_mxy, e_my), tspan_meta,frac_occup_init,options2);
    %frac_occup(end,1) = 0;
    
    occupancy_del_m(i, :) = frac_occup(end,:);
    cmx_collector(i,:) = c_mx;
    
    
end
%save ("new_params_coexist_nonneg_allcombos.mat")
%% Figures
figure()
plot(del_m, occupancy_del_m(:, 1:3))
xlabel('mutualist dispersal rate')
%xlim([1.0 29.0]);
ylabel('fraction of patches occupied')
%xline([1.0 2.6, 10.4, 26.9],'--',{'Exploitative (x extinct)','Mutualism (y fitter)','Mutualism (x fitter)', 'Mutualism (y fitter)'})
title('Fraction of patches occupied vs mutualist dispersal rate')
legend('Species with mutualist (x)', 'Species without mutualist (y)', 'mutualist (m)', 'location', 'best' )
%fig1name = sprintf('occupancy_vs_del_m.jpeg');
%print('occupancy_vs_del_m_allc','-djpeg','-r600')






