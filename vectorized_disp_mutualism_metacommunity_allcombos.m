%%%%%%%% This code vectorizes the ODEs that need to be solved for multiple
%%%%%%%% varying parameters (here 'q' and 'del_m')
clear all
%% parameter definitions
r_x = 5; r_y = 5;

alpha_xy = 0.73; alpha_yx = 0.60; %really only affects y's local density and dispersal

K_x = 200; K_y = 200; %no point touching this

del_x = 0.01; del_y = 0.05;
%q = 0.45;
d_m = 0.01; a = 0.8;%reducing 'a' reduces dispersal rates where px > py
del_m = 3:1:53; %I'm only not starting from zero because the computational costs are absurd
%a = 0.51:0.01:0.81;%0.7-1.1 seems to work for this fig.
q = 0.73:0.01:0.80;
[Del,Q] = meshgrid(del_m,q);

Del_col = Del(:); Q_col = Q(:); 

num_combinations = numel(Q_col); 

k_eff = 1.0; %efficiency of dispersing seeds to habitable patches

z_x = 0.7; z_y = 0.7; z_m = 0.4; %scaling factors for patch extinction rates. changing z_m relative to z_x and z_m does not change qual. change results

e_xmin = 0.05;e_ymin = 0.05; e_mmin = 0.3; e_mxmin = e_xmin; 

tspan = [0,1000];

k_x = 0.08; k_y = 0.08; k_m = 0.10; %

x_init = 0.1; y_init = 0.1; m_init = 0.1;

spp_init_no_m = [x_init; y_init; 0];
spp_init_no_y = [x_init; 0; m_init];
spp_init_no_x = [0; y_init; m_init];
spp_init = [x_init; y_init; m_init];

f = 1; % f is the fraction of diet consumed by frugivore that consists of 'x'. (1-f) is fraction that is 'y'
%variable collectors across parameter sweeps

% occupancy_del_m = zeros(size(del_m,1),3);
% eta_check = zeros(size(del_m,1),1);
% mu_collector = zeros(size(del_m,1),1);
% gamma = zeros(size(del_m,1),1);
% cm_collector = zeros(size(del_m,1),1);
% eta2 = zeros(size(del_m,1),1);
% em_collector = zeros(size(del_m,1),1);
% lambda_collector = zeros(size(del_m,1),1);

%% Local patch dynamics
threshold = 10^-8;
options = odeset('NonNegative',[1,2,3]);%,'Events',@nonNegativeEvent);
% [t_patch_no_m,local_dens_no_m] = ode45(@(t,y) vectorized_LocalSpeciesInteraction(t,y,r_x,r_y,alpha_xy,alpha_yx, K_x, K_y, del_x, del_y, Del_col', a, Q_col', d_m, num_combinations), tspan./10, repmat(spp_init_no_m,1,num_combinations));
% [t_patch_no_y,local_dens_no_y] = ode45(@(t,y) vectorized_LocalSpeciesInteraction(t,y,r_x,r_y,alpha_xy,alpha_yx, K_x, K_y, del_x, del_y, Del_col', a, Q_col', d_m, num_combinations), tspan./10, repmat(spp_init_no_y,1,num_combinations));
% [t_patch,local_dens] = ode45(@(t,y)vectorized_LocalSpeciesInteraction(t,y,r_x,r_y,alpha_xy,alpha_yx, K_x, K_y, del_x, del_y,Del_col',a, Q_col', d_m, num_combinations), tspan, repmat(spp_init, 1, num_combinations), options);

% in case the frugivore is a generalist (or specialist if f=1)
[t_patch_no_m,local_dens_no_m] = ode45(@(t,y) vectorized_LocalSpeciesInteraction_generalist(t,y,r_x,r_y,alpha_xy,alpha_yx, K_x, K_y, del_x, del_y, Del_col', a, Q_col', d_m, num_combinations, f), tspan, repmat(spp_init_no_m,1,num_combinations),options);
[t_patch_no_y,local_dens_no_y] = ode45(@(t,y) vectorized_LocalSpeciesInteraction_generalist(t,y,r_x,r_y,alpha_xy,alpha_yx, K_x, K_y, del_x, del_y, Del_col', a, Q_col', d_m, num_combinations, f), tspan, repmat(spp_init_no_y,1,num_combinations), options);
[t_patch,local_dens] = ode45(@(t,y)vectorized_LocalSpeciesInteraction_generalist(t,y,r_x,r_y,alpha_xy,alpha_yx, K_x, K_y, del_x, del_y,Del_col',a, Q_col', d_m, num_combinations, f), tspan, repmat(spp_init, 1, num_combinations), options);
% while t_patch(end) < tspan(end)
%     % Update initial conditions based on the last state
%     ye(abs(ye(3:3:end)) <= threshold) = 0;
%     initialConditions = ye'
%     % Solve again starting from the last time point
%     [t_temp,y_temp,te,ye] = ode45(@(t,y) vectorized_LocalSpeciesInteraction_generalist(t,y,r_x,r_y,alpha_xy,alpha_yx, K_x, K_y, del_x, del_y, Del_col', a, Q_col', d_m, num_combinations, f), [te tspan(end)], initialConditions, options);
%     % Concatenate the results
%     t_patch = [t_patch; t_temp(2:end)];  % Avoid duplicating the common time point
%     local_dens = [local_dens; y_temp(2:end, :)];
% end

num_timepoints_no_m = length(t_patch_no_m);
num_timepoints_no_y = length(t_patch_no_y);
num_timepoints = length(t_patch);
 
local_dens_no_m_3d = reshape(local_dens_no_m(end,:),[],num_combinations);

local_dens_no_y_3d = reshape(local_dens_no_y(end,:),[],num_combinations);
local_dens_no_y_3d(3, local_dens_no_y_3d (3,:) < threshold) = 0; %to prevent torture by numerically minuscule values that are biologically irrelevant

local_dens_3d = reshape(local_dens(end,:),[],num_combinations);
local_dens_3d(3, local_dens_3d (3,:) < threshold) = 0;




%      figure()
%      plot(t_patch,local_dens)
%% metacommunity dynamics, patch colonisation/extinction rates
% e_x = e_xmin..*((K_x../(local_dens_no_m_3d(1,:)))..^z_x);e_y = e_ymin..*(K_y../(local_dens_3d(2,:)))..^z_y;
% e_mx = e_xmin..*(K_x../(local_dens_3d(1,:)))..^z_x;
% e_m = e_mmin..*(K_x../(local_dens_3d(3,:)))..^z_m; %assume max population size of mutualist is similar to that of x and y

e_x0 = e_xmin;e_y0 = e_ymin;
e_xy = e_xmin.*((K_x./(local_dens_no_m_3d(1,:))).^z_x); e_yx = e_ymin.*((K_y./(local_dens_no_m_3d(2,:))).^z_y);
e_xm = e_xmin.*((K_x./(local_dens_no_y_3d(1,:))).^z_x); e_ym = e_ymin;
e_xym = e_xmin.*(K_x./(local_dens_3d(1,:))).^z_x; e_yxm = e_ymin.*(K_y./(local_dens_3d(2,:))).^z_y;
e_mx = e_mmin.*(K_x./(local_dens_no_y_3d(3,:))).^z_m; e_my = 0; e_mxy = e_mmin.*(K_x./(local_dens_3d(3,:))).^z_m; %assume max population size of mutualist is similar to that of x and y
    
%mu = e_mx - e_x;

tspan_meta = [0,1000];

% c_x = (k_x.*del_x)..*local_dens_no_m_3d(1,:);
% c_y = k_y.*del_y..*local_dens_3d(2,:);
% c_m = k_m.*Del_col'..*local_dens_3d(3,:);%repelem(del_m,numel(q)) can also be used in place of Del_col' if u feel funky
% c_mx =(k_x.*del_x+k_m.*k_eff.*a..*(1-Q_col')..*Del_col'..*local_dens_3d(3,:))..*local_dens_3d(1,:);

%c_x0 = (k_x.*del_x).*K_x; c_y0 = (k_y.*del_y).*K_y; %patches with only one species
c_x0 = (k_x.*del_x).*K_x*(1-(del_x/r_x));
c_y0 = (k_y.*del_y).*K_y*(1-(del_y/r_y)); 
c_xy = (k_x.*del_x).*local_dens_no_m_3d(1,:); c_yx = (k_y.*del_y).*local_dens_no_m_3d(2,:); c_mx = k_m.*Del_col'.*local_dens_no_y_3d(3,:); c_my = 0;%patches with 2 species
c_xm = ((k_x.*del_x)+k_m*k_eff.*Del_col'.*local_dens_no_y_3d(3,:)).*local_dens_no_y_3d(1,:);%k_x.*del_x.*local_dens_no_y_3d(1,:); 
c_ym = (k_y.*del_y).*K_y*(1-(del_y/r_y));%(k_y.*del_y).*K_y; % patches with one plant-one frugivore
c_xym = (k_x.*del_x+k_m.*k_eff.*Del_col'.*local_dens_3d(3,:)).*local_dens_3d(1,:); c_yxm = k_y.*del_y.*local_dens_3d(2,:); c_mxy = k_m.*Del_col'.*local_dens_3d(3,:);%all species present


%lambda = c_mx-c_x;
%% metacommunity equations
frac_occup_init = [0.1,0.1,0.1];
IC_betweenpatch = repelem(frac_occup_init,num_combinations);
IC_betweenpatch = reshape(IC_betweenpatch,[],num_combinations);
options2 = odeset('NonNegative',[1,2,3]);
if any(e_mxy == inf) || any(e_mx ==inf)
    IC_betweenpatch(3,find(e_mx == inf)) = 0;
    IC_betweenpatch(3,find(e_mxy == inf)) = 0;
end

[t_syst, frac_occup] = ode45(@(t,y)vectorized_BetweenPatchDynamics_allcombos(t,y, c_x0, c_xm, c_xym, c_xy, c_y0, c_ym, c_yxm, c_yx, c_mx, c_mxy, c_my, e_x0, e_xy, e_xm, e_xym, e_y0, e_yx, e_ym, e_yxm, ...
    e_mx, e_mxy, e_my, num_combinations), tspan_meta, IC_betweenpatch(:), options2);


num_tpts_t_syst = length(t_syst); frac_occup_3d = reshape(frac_occup(end,:), [],num_combinations);
% 
% eta_check = 0.5..*sqrt((((lambda..*(1+(e_m../c_m))+c_x)-(e_x+mu))../(c_x+lambda))..^2 + 4..*(e_m../c_m)..*((mu-lambda)../(c_x+lambda)));
%      figure()
%      plot(t_syst(end-100:end),frac_occup(end-100:end,1:2));

% mu_collector(i,1) = mu;
% cm_collector(i,1) = c_m;
% em_collector(i,1) = e_m;
% lambda_collector(i,1) = lambda;
% occupancy_del_m(i, :) = frac_occup(end,:);
% gamma(i,1) = (0.5..*((lambda.*(1+(e_m../c_m))+c_x)-(e_x+mu)))../(c_x+lambda);%0.5.*(1-((mu-e_x-e_m)./(c_x+c_m)));
% eta2(i,:) = 4.*(e_m../c_m).*((mu-lambda)../(c_x+lambda));
%save ("vectorized_q_del_m_varied.mat")
%save(strcat(['occupancy_qlow_' num2str(q(1)*100, '%d') '_qhi_' num2str(q(end)*100, '%d') '_delmlo_' num2str(del_m(1), '%d') '_delmhi_' num2str(del_m(end), '%d')]), Del,Q,frac_occup_3d,q,del_m, '-mat');

save (sprintf('occupancy_qlow_%d_qhi_%d_delmlo_%d_delmhi_%d.mat',q(1)*100,q(end)*100,del_m(1),del_m(end)),'Del','Q','frac_occup_3d','q','del_m');

% save(filename,Del,Q,frac_occup_3d,q,del_m)
%% surface plot
figure()
surf(Del, Q, reshape(frac_occup_3d(1,:),numel(q),numel(del_m)))
colormap parula; freezeColors;
hold on
surf(Del, Q, reshape(frac_occup_3d(2,:),numel(q),numel(del_m)))
surf(Del, Q, reshape(frac_occup_3d(3,:),numel(q),numel(del_m)))
colormap winter; freezeColors;
%surf(Del, Q,occupancy_del_m(:,:,1))
xlabel('mutualist dispersal rate (\delta_m)')
ylabel ('consumption fraction (q)')
%xlim([1.0 29.0]);
zlabel('fraction of patches occupied')
%xline([1.0 2.6, 10.4, 26.9],'--',{'Exploitative (x extinct)','Mutualism (y fitter)','Mutualism (x fitter)', 'Mutualism (y fitter)'})
title('Fraction of patches occupied vs mutualist dispersal rate and predation rate')
legend('Species with mutualist (x)')%, 'Species without mutualist (y)', 'mutualist (m)', 'location', 'best' )
%fig1name = sprintf('occupancy_vs_del_m.jpeg');
print('specialist_occup_q_vs_del_m','-djpeg','-r600')
