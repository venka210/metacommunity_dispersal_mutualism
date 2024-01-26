clear vars

%% parameter definition
r_x = 5; r_y = 5;

alpha_xy = 0.73; alpha_yx = 0.60; %really only affects y's local density and dispersal

K_x = 200; K_y = 200; %no point touching this

del_x = 0.01; del_y = 0.05; del_m = 3:1:27; %I'm only not starting from zero because the computational costs are absurd

a = 0.70:0.01:0.90;

q = 0.45; d_m = 0.03; %

k_eff = 1; %efficiency of dispersing seeds to habitable patches

z_x = 0.7; z_y = 0.7; z_m = 0.4; %scaling factors for patch extinction rates. changing z_m relative to z_x and z_m does not change qual. change results

e_xmin = 0.05;e_ymin = 0.05; e_mmin = 0.05; e_mxmin = e_xmin;

tspan = [0,1000];

k_x = 0.08; k_y = 0.08; k_m = 0.08;

x_init = 0.1; y_init = 0.1; m_init = 0.1;
    
spp_init_no_m = [x_init; y_init; 0];
spp_init = [x_init; y_init; m_init];

%variable collectors across parameter sweeps

occupancy_del_m = zeros(size(del_m,1),size(a,1),3);
eta_check = zeros(size(del_m,1),size(a,1),1);
mu_collector = zeros(size(del_m,1),size(a,1),1);
gamma = zeros(size(del_m,1),size(a,1),1);
cm_collector = zeros(size(del_m,1),size(a,1),1);
eta2 = zeros(size(del_m,1),size(a,1),1);
em_collector = zeros(size(del_m,1),size(a,1),1);
lambda_collector = zeros(size(del_m,1),size(a,1),1);


%% within patch dynamics
for i = 1:length(del_m)
    for j = 1:length(a)
        
        [t_patch_no_m,local_dens_no_m] = ode45(@(t,y)LocalSpeciesInteraction(t,y,r_x,r_y,alpha_xy,alpha_yx, K_x, K_y, del_x, del_y, del_m(i), a(j), q, d_m), tspan/10, spp_init_no_m);
        [t_patch,local_dens] = ode45(@(t,y)LocalSpeciesInteraction(t,y,r_x,r_y,alpha_xy,alpha_yx, K_x, K_y, del_x, del_y, del_m(i), a(j), q, d_m), tspan, spp_init);

    %      figure()
    %      plot(t_patch,local_dens)
        %% metacommunity dynamics

        e_x = e_xmin*((K_x/(local_dens_no_m(end,1)))^z_x);e_y = e_ymin*(K_y/(local_dens(end,2)))^z_y;
        e_mx = e_xmin*(K_x/(local_dens(end,1)))^z_x; e_mxy=0;
        e_m = e_mmin*(K_x/(local_dens(end,3)))^z_m; %assume max population size of mutualist is similar to that of x and y

    %     if local_dens(end,3) == 0
    %         e_mx = 0;
    %     else
    %         e_mx = e_xmin*(K_x^z_m)*(local_dens(end,3)^(-z_m)); %density dependent patch extinction rates
    %     end

    %     if e_m == 0
    %         mu = 0;
    %     else
    %         mu = e_mx - e_x;
    %     end
        mu = e_mx - e_x;

        tspan_meta = [0, 1000];

        c_x = (k_x*del_x)*local_dens_no_m(end,1);
        c_y = k_y*del_y*local_dens(end,2);
        c_m = k_m*del_m(i)*local_dens(end,3);
        c_mx =(k_x*del_x+k_m*k_eff*a(j)*(1-q)*del_m(i)*local_dens(end,3))*local_dens(end,1); c_mxy = 0;

        lambda = c_mx-c_x;
        options = odeset('NonNegative',1);

        [t_syst, frac_occup] = ode15s(@(t,y)BetweenPatchDynamics_new(t,y, c_x, c_y, c_m, e_x, e_y, e_m, mu, lambda), tspan_meta, (local_dens(end,:))',options);
    %     [t_syst, frac_occup] = ode15s(@(t,y)BetweenPatchDynamics_allcombos(t,y, c_x0, c_xm, c_xym, c_xy, c_y0, c_ym, c_yxm, c_yx, c_mx, c_mxy, e_x0, e_xy, e_xm, e_xym, e_y0, e_yx, e_ym, e_yxm, ...
    % e_mx, e_mxy);
        %frac_occup(end,1) = 0;


        %eta_check(i,1) = (e_m/c_m)*((c_m + mu)/(c_m+c_x));
        eta_check(i,j,1) = 0.5*sqrt((((lambda*(1+(e_m/c_m))+c_x)-(e_x+mu))/(c_x+lambda))^2 + 4*(e_m/c_m)*((mu-lambda)/(c_x+lambda)));
        %0.5*sqrt(((1-((mu-e_x-e_m)/(c_x+c_m)))^2)+(4*(e_m/c_m)*((mu-c_m)/(c_x+c_m)))); %for bascompte model to prevent imaginary roots
    %      figure()
    %      plot(t_syst(end-100:end),frac_occup(end-100:end,1:2));

        mu_collector(i,j,1) = mu;
        cm_collector(i,j,1) = c_m;
        em_collector(i,j,1) = e_m;
        lambda_collector(i,j,1) = lambda;
        occupancy_del_m(i,j,1) = frac_occup(end,1);
        occupancy_del_m(i,j,2) = frac_occup(end,2);
        occupancy_del_m(i,j,3) = frac_occup(end,3);
        gamma(i,j,1) = (0.5*((lambda*(1+(e_m/c_m))+c_x)-(e_x+mu)))/(c_x+lambda);%0.5*(1-((mu-e_x-e_m)/(c_x+c_m)));
        eta2(i,j,:) = 4*(e_m/c_m)*((mu-lambda)/(c_x+lambda));

    %     del_m(i)
    %     if mu > lambda
    %         disp('positive root only allowed')
    %         disp(frac_occup(end,1))
    %     else
    %         if (gamma(i,1)^2) > eta2(i,:)
    %             disp('both roots allowed')
    %             disp(frac_occup(end,1))
    %         else
    %             disp('no roots allowed')
    %         end
    %     end


    end
end
%save ("new_params_coexist_nonneg.mat")
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
print('occupancy_vs_del_m','-djpeg','-r600')
%%
% figure()
% plot(del_m, gamma+eta_check)
% hold on
% plot(del_m, occupancy_del_m(:,1),'o')
% %plot(del_m, occupancy_del_m(:,3),'o')
% plot(del_m, gamma-eta_check)
% hold off
% xlabel('mutualist dispersal rate (\delta_m)')
% ylabel('Roots of {p_x}^*')
% %xlim([1.0 29.0]);
% %ylim([-1.05 1.05]);
% %xline([1.0 2.6, 10.4, 26.9],'--',{'Exploitative (x extinct)','Mutualism (y fitter)','Mutualism (x fitter)', 'Mutualism (y fitter)'})
% %yline(0.5*(1-((mu-e_x-e_m)/(c_x+c_m))));
% %yline(0.00)
% fig1name = sprintf('eta_vs_del_m.jpeg');
% print(fig1name,'-djpeg','-r600')
%%
figure()
[Del,A] = meshgrid(del_m,a);
surf(Del', A', occupancy_del_m(:,:,1))
%surf(Del, A,occupancy_del_m(:,:,1))
xlabel('mutualist dispersal rate')
ylabel ('predation rate')
%xlim([1.0 29.0]);
zlabel('fraction of patches occupied')
%xline([1.0 2.6, 10.4, 26.9],'--',{'Exploitative (x extinct)','Mutualism (y fitter)','Mutualism (x fitter)', 'Mutualism (y fitter)'})
title('Fraction of patches occupied vs mutualist dispersal rate and predation rate')
legend('Species with mutualist (x)')%, 'Species without mutualist (y)', 'mutualist (m)', 'location', 'best' )
%fig1name = sprintf('occupancy_vs_del_m.jpeg');
print('occupancy_3d_vs_del_m','-djpeg','-r600')







