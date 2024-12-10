clear all

%% parameter definition
r_x = 5; r_y = 5;

alpha_xy = 0.73; alpha_yx = 0.60; %really only affects y's local density and dispersal -- alpha_xy = 0.73; alpha_yx = 0.60;

K_x = 200; K_y = 200; %no point touching this

del_x = 0.01; del_y = 0.03; del_m = 0:0.1:10;%del_m = 0:0.1:4; %I'm only not starting from zero because the computational costs are absurd

k_eff = 0.1; %efficiency of dispersing seeds to habitable patches

z_x = 0.9; z_y = 0.9; z_m = 0.6; %scaling factors for patch extinction rates. changing z_m relative to z_x and z_m does not change qual. change results

e_xmin = 0.05;e_ymin = 0.05; e_mmin = 0.008; %e_mxmin = 0.01;

tspan = [0,200];

k_x = 0.1; k_y = 0.1; k_m = 0.03;

d_m = 1; a = 1.0; %reducing 'a' reduces dispersal rates where px > py
q = 1.0; 
x_init = 0.1; y_init = 0.1; m_init = 0.1;
    
spp_init_no_m = [x_init; y_init; 0];
spp_init_no_y = [x_init; 0; m_init];
spp_init_no_x = [0; y_init; m_init];
spp_init = [x_init; y_init; m_init];

%variable collectors across parameter sweeps

occupancy_del_m = zeros(size(del_m,1),3);



%% within patch dynamics

for i = 1:length(del_m)
    
    options = odeset('NonNegative',[1,2,3]);
    [t_patch_no_m,local_dens_no_m] = ode45(@(t,y)LocalSpeciesInteraction(t,y,r_x,r_y,alpha_xy,alpha_yx, K_x, K_y, del_x, del_y, del_m(i), a, q, d_m), tspan/10, spp_init_no_m, options);
    [t_patch_no_y,local_dens_no_y] = ode45(@(t,y)LocalSpeciesInteraction(t,y,r_x,r_y,alpha_xy,alpha_yx, K_x, K_y, del_x, del_y, del_m(i), a, q, d_m), tspan/10, spp_init_no_y,options);
    [t_patch,local_dens] = ode45(@(t,y)LocalSpeciesInteraction(t,y,r_x,r_y,alpha_xy,alpha_yx, K_x, K_y, del_x, del_y, del_m(i), a, q, d_m), tspan, spp_init, options);

    %      figure()
    %      plot(t_patch,local_dens)
        %% metacommunity dynamics
    e_x0 = e_xmin;e_y0 = e_ymin;
    e_xy = e_xmin.*((K_x./(local_dens_no_m(end,1))).^z_x); e_yx = e_ymin.*((K_y./(local_dens_no_m(end,2))).^z_y);
    e_xm = e_xmin.*((K_x./(local_dens_no_y(end,1))).^z_x); e_ym = e_ymin;
    e_mx = e_mmin.*(K_x./(local_dens_no_y(end,3))).^z_m;% remove to checkno mutualism case
    e_my = 0;
    %e_mx = 0; %in case we want to analyse no mutualism setting
    %e_mxy = 0;%ic no mutualism
    e_mxy = e_mmin.*(K_x./(local_dens(end,3))).^z_m; 
    e_xym = e_xmin.*(K_x./(local_dens(end,1))).^z_x; e_yxm = e_ymin.*(K_y./(local_dens(end,2))).^z_y;%assume max population size of mutualist is similar to that of x and y
    
    tspan_meta = [0,500];
    
    %patches with only one species
    c_x0 = (k_x.*del_x).*K_x*(1-(del_x/r_x));
    c_y0 = (k_y.*del_y).*K_y*(1-(del_y/r_y)); 
    
    %patches with 2 species
    c_xy = (k_x.*del_x).*local_dens_no_m(end,1); 
    c_xm = (k_x*del_x).*local_dens_no_y(end,1)+k_m*k_eff.*del_m(i).*local_dens_no_y(end,3).*(K_x - local_dens_no_y(end,1));
    
    c_yx =(k_y.*del_y).*local_dens_no_m(end,2);
    c_ym = (k_y.*del_y).*K_y*(1-(del_y/r_y)); 
    
    c_mx = k_m.*del_m(i).*local_dens_no_y(end,3); 
    c_my = 0;
    
    %all species present
    c_xym = k_x*del_x.*local_dens(end,1)+k_m*k_eff.*del_m(i).*local_dens(end,3).*(K_x - local_dens(end,1));
    c_yxm = k_y*del_y*local_dens(end,2);
    c_mxy = k_m.*del_m(i).*local_dens(end,3);%all species present

    
    
    %substitute local_dens_3d(*spp ID*,:) for local_dens(end, *spp ID*) in non-vectorized code
    %lambda = c_mx-c_x;
    options2 = odeset('NonNegative',[1,2,3]);
    frac_occup_init = [0.01;0.01;0.01]; %initial patch occupancies of each species x, y, and, m respectively. 
    % if e_mxy == inf
    %     frac_occup_init(3,1) = 0;
    % end

    [t_syst, frac_occup] = ode45(@(t,y)vectorized_BetweenPatchDynamics_allcombos(t,y, c_x0, c_xm, c_xym, c_xy, c_y0, c_ym, c_yxm, c_yx, c_mx, c_mxy, c_my, e_x0, e_xy, e_xm, e_xym, e_y0, e_yx, e_ym, e_yxm, ...
    e_mx, e_mxy, e_my, 1), tspan_meta, frac_occup_init, options2);
    %frac_occup(end,1) = 0;
    
    
    %eta_check(i,1) = (e_m/c_m)*((c_m + mu)/(c_m+c_x));
%     eta_check(i,1) = 0.5*sqrt((((lambda*(1+(e_m/c_m))+c_x)-(e_x+mu))/(c_x+lambda))^2 + 4*(e_m/c_m)*((mu-lambda)/(c_x+lambda)));
%     %0.5*sqrt(((1-((mu-e_x-e_m)/(c_x+c_m)))^2)+(4*(e_m/c_m)*((mu-c_m)/(c_x+c_m)))); %for bascompte model to prevent imaginary roots
% %      figure()
% %      plot(t_syst(end-100:end),frac_occup(end-100:end,1:2));
% 
%     mu_collector(i,1) = mu;
%     cm_collector(i,1) = c_m;
%     em_collector(i,1) = e_m;
%     lambda_collector(i,1) = lambda;
     occupancy_del_m(i, :) = frac_occup(end,:);

    
end
save (sprintf('specialist_occupancy_q_%0.2f_del_m_max_%3.1f.mat',q,del_m(end)));
%save (sprintf('specialist_occupancy_no_m_q_%0.2f_del_m_max_%3.1f.mat',q,del_m(end)));
%% Figures

figure()
plot(t_syst,frac_occup(:,1),'LineWidth',2); hold on; plot(t_syst, frac_occup(:,2),'LineWidth',2); plot(t_syst, frac_occup(:,3),'--','LineWidth',2)
xlabel('time (global dynamics)')
ylabel('fraction of patches occupied')
title('Fraction of patches occupied vs time')
legend('Species with mutualist (x)', 'Species without mutualist (y)', 'mutualist (m)', 'location', 'best' )
print(sprintf('occupancy_vs_time_q_%0.2f_del_m_max_%3.1f.jpg',q,del_m(end)),'-djpeg','-r600')
%print(sprintf('occupancy_vs_time_no_m_q_%0.2f_del_m_max_%3.1f.jpg',q,del_m(end)),'-djpeg','-r600')

figure()
plot(t_patch,local_dens(:,1),'LineWidth',2); hold on; plot(t_patch, local_dens(:,2),'LineWidth',2); plot(t_patch, local_dens(:,3),'--','LineWidth',2)
xlabel('time (local dynamics)')
ylabel('local patch density')
title('Local patch population density vs time')
legend('Species with mutualist (x)', 'Species without mutualist (y)', 'mutualist (m)', 'location', 'best' )
print(sprintf('locdens_vs_time_q_%0.2f_del_m_max_%3.1f.jpg',q,del_m(end)),'-djpeg','-r600')
%print(sprintf('locdens_vs_time_no_m_q_%0.2f_del_m_max_%3.1f.jpg',q,del_m(end)),'-djpeg','-r600')

figure()
plot(del_m, occupancy_del_m(:, 1:3),'LineWidth',3)
xlabel('mutualist dispersal rate (\delta_m)')
yline(0.1368,'g--','LineWidth',3)
%xlim([1.0 29.0]);
ylabel('fraction of patches occupied (p_x, p_y, p_m)')
%xline([1.0 2.6, 10.4, 26.9],'--',{'Exploitative (x extinct)','Mutualism (y fitter)','Mutualism (x fitter)', 'Mutualism (y fitter)'})
title('Fraction of patches occupied vs mutualist dispersal rate')
legend('Species with mutualist (x)', 'Species without mutualist (y)', 'mutualist (m)', 'p_x w/o mutualist','location', 'east')
%fig1name = sprintf('occupancy_vs_del_m.jpeg');
print(sprintf('occupancy_vs_delm_q_%0.2f_del_m_max_%3.1f.jpg',q,del_m(end)),'-djpeg','-r600')
%print(sprintf('occupancy_vs_delm_no_m_q_%0.2f_del_m_max_%3.1f.jpg',q,del_m(end)),'-djpeg','-r600')







