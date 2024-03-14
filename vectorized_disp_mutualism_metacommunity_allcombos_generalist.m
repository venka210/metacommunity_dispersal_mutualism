
%%%%%%%% This code vectorizes the ODEs that need to be solved for multiple
%%%%%%%% varying parameters (here 'q' and 'del_m')
clear all
%% parameter definitions

r_x = 5; r_y = 5;

alpha_xy = 0.73; alpha_yx = 0.60; %really only affects y's local density and dispersal -- alpha_xy = 0.73; alpha_yx = 0.60;

K_x = 200; K_y = 200; %no point touching this

del_x = 0.01; del_y = 0.03; %del_m = 0:0.1:4; %I'm only not starting from zero because the computational costs are absurd

%a = 1; q = 1; d_m = 1;

k_eff = 0.1; %efficiency of dispersing seeds to habitable patches

z_x = 0.8; z_y = 0.8; z_m = 0.5; %scaling factors for patch extinction rates. changing z_m relative to z_x and z_m does not change qual. change results

e_xmin = 0.05;e_ymin = 0.05; e_mmin = 0.05; e_mxmin = e_xmin;

tspan = [0,200];

k_x = 0.1; k_y = 0.1; k_m = 0.05;

%x_init = 0.1; y_init = 0.1; m_init = 0.1;

%del_x = 0.01; del_y = 0.05;

%q = 0.45;
d_m = 1; a = 1.0;%reducing 'a' reduces dispersal rates where px > py

del_m_all = 0:0.1:12; 

q_all = 0.0:0.01:1.0;

all_params = [repelem(del_m_all,size(q_all,2));repmat(q_all,1,size(del_m_all,2))];

f = 0.96; % f is the fraction of diet consumed by frugivore that consists of 'x'. (1-f) is fraction that is 'y'

[Del_and_Q, invalid_ind, good_bad] = parameterRangeCheck(del_m_all,q_all,f);

%del_m = unique(Del_and_Q(1,:)); q = unique(Del_and_Q(2,:));

%[allDel,allQ] = meshgrid(del_m_all,q_all);

Del = Del_and_Q(1,:); Q = Del_and_Q(2,:);

Del_col = Del(:); Q_col = Q(:); 

num_combinations = size(Del_and_Q,2);   %numel(Q_col); 
%%
%k_eff = 1; %efficiency of dispersing seeds to habitable patches

%z_x = 0.7; z_y = 0.7; z_m = 0.5; %scaling factors for patch extinction rates. changing z_m relative to z_x and z_m does not change qual. change results

%e_xmin = 0.05;e_ymin = 0.05; e_mmin = 0.03; e_mxmin = e_xmin; 

%tspan = [0,1000];

%k_x = 0.08; k_y = 0.08; k_m = 0.02; %

x_init = 0.01; y_init = 0.01; m_init = 0.01;

spp_init_no_m = [x_init; y_init; 0];
spp_init_no_y = [x_init; 0; m_init];
spp_init_no_x = [0; y_init; m_init];
spp_init = [x_init; y_init; m_init];

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
thresholdValue = 10^-6;
%options = odeset('Events',@nonNegativeEvent); %options2 =  odeset('Events',@(t, y) eventFunction(t, y, thresholdValue));
options = odeset('NonNegative',[1,2,3]);
% [t_patch_no_m,local_dens_no_m] = ode45(@(t,y) vectorized_LocalSpeciesInteraction(t,y,r_x,r_y,alpha_xy,alpha_yx, K_x, K_y, del_x, del_y, Del_col', a, Q_col', d_m, num_combinations), tspan./10, repmat(spp_init_no_m,1,num_combinations));
% [t_patch_no_y,local_dens_no_y] = ode45(@(t,y) vectorized_LocalSpeciesInteraction(t,y,r_x,r_y,alpha_xy,alpha_yx, K_x, K_y, del_x, del_y, Del_col', a, Q_col', d_m, num_combinations), tspan./10, repmat(spp_init_no_y,1,num_combinations));
% [t_patch,local_dens] = ode45(@(t,y)vectorized_LocalSpeciesInteraction(t,y,r_x,r_y,alpha_xy,alpha_yx, K_x, K_y, del_x, del_y,Del_col',a, Q_col', d_m, num_combinations), tspan, repmat(spp_init, 1, num_combinations), options);

% in case the frugivore is a generalist
[t_patch_no_x,local_dens_no_x] = ode45(@(t,y) vectorized_LocalSpeciesInteraction_generalist(t,y,r_x,r_y,alpha_xy,alpha_yx, K_x, K_y, del_x, del_y, Del_col', a, Q_col', d_m, num_combinations, f), tspan, repmat(spp_init_no_x,1,num_combinations), options);
[t_patch_no_m,local_dens_no_m] = ode45(@(t,y) vectorized_LocalSpeciesInteraction_generalist(t,y,r_x,r_y,alpha_xy,alpha_yx, K_x, K_y, del_x, del_y, Del_col', a, Q_col', d_m, num_combinations, f), tspan, repmat(spp_init_no_m,1,num_combinations), options);
[t_patch_no_y,local_dens_no_y] = ode45(@(t,y) vectorized_LocalSpeciesInteraction_generalist(t,y,r_x,r_y,alpha_xy,alpha_yx, K_x, K_y, del_x, del_y, Del_col', a, Q_col', d_m, num_combinations, f), tspan, repmat(spp_init_no_y,1,num_combinations), options);
[t_patch,local_dens] = ode45(@(t,y)vectorized_LocalSpeciesInteraction_generalist(t,y,r_x,r_y,alpha_xy,alpha_yx, K_x, K_y, del_x, del_y,Del_col',a, Q_col', d_m, num_combinations, f), tspan, repmat(spp_init, 1, num_combinations), options);

num_timepoints_no_x = length(t_patch_no_x);
num_timepoints_no_m = length(t_patch_no_m);
num_timepoints_no_y = length(t_patch_no_y);
num_timepoints = length(t_patch);

local_dens_no_x_3d = reshape(local_dens_no_x(end,:),[],num_combinations);
local_dens_no_x_3d(3, local_dens_no_x_3d (3,:) < thresholdValue) = 0;

local_dens_no_m_3d = reshape(local_dens_no_m(end,:),[],num_combinations);

local_dens_no_y_3d = reshape(local_dens_no_y(end,:),[],num_combinations);
local_dens_no_y_3d(3, local_dens_no_y_3d (3,:) < thresholdValue) = 0; %to prevent torture by numerically minuscule values that are biologically irrelevant

local_dens_3d = reshape(local_dens(end,:),[],num_combinations);
local_dens_3d(3, local_dens_3d (3,:) < thresholdValue) = 0;

%      figure()
%      plot(t_patch,local_dens)
%% metacommunity dynamics
% e_x = e_xmin..*((K_x../(local_dens_no_m_3d(1,:)))..^z_x);e_y = e_ymin..*(K_y../(local_dens_3d(2,:)))..^z_y;
% e_mx = e_xmin..*(K_x../(local_dens_3d(1,:)))..^z_x;
% e_m = e_mmin..*(K_x../(local_dens_3d(3,:)))..^z_m; %assume max population size of mutualist is similar to that of x and y

e_x0 = e_xmin;e_y0 = e_ymin;
e_xy = e_xmin.*((K_x./(local_dens_no_m_3d(1,:))).^z_x); e_yx = e_ymin.*((K_y./(local_dens_no_m_3d(2,:))).^z_y);
e_xm = e_xmin.*((K_x./(local_dens_no_y_3d(1,:))).^z_x); e_ym = e_ymin.*((K_y./(local_dens_no_x_3d(2,:))).^z_y);
e_xym = e_xmin.*(K_x./(local_dens_3d(1,:))).^z_x; e_yxm = e_ymin.*(K_y./(local_dens_3d(2,:))).^z_y;
e_mx = e_mmin.*(K_x./(local_dens_no_y_3d(3,:))).^z_m; e_my = e_mmin.*(K_x./(local_dens_no_x_3d(3,:))).^z_m; 
e_mxy = e_mmin.*(K_x./(local_dens_3d(3,:))).^z_m; %assume max population size of mutualist is similar to that of x and y
    
%mu = e_mx - e_x;

tspan_meta = [0,1000];

% c_x = (k_x.*del_x)..*local_dens_no_m_3d(1,:);
% c_y = k_y.*del_y..*local_dens_3d(2,:);
% c_m = k_m.*Del_col'..*local_dens_3d(3,:);%repelem(del_m,numel(q)) can also be used in place of Del_col' if u feel funky
% c_mx =(k_x.*del_x+k_m.*k_eff.*a..*(1-Q_col')..*Del_col'..*local_dens_3d(3,:))..*local_dens_3d(1,:);

%patches with only one species
c_x0 = (k_x.*del_x).*K_x*(1-(del_x/r_x));
c_y0 = (k_y.*del_y).*K_y*(1-(del_y/r_y)); 

%patches with 2 species
c_xy = (k_x.*del_x).*local_dens_no_m_3d(1,:); 
c_xm = (k_x*del_x).*local_dens_no_y_3d(1,:)+k_m*k_eff.*Del_col'.*local_dens_no_y_3d(3,:).*(K_x - local_dens_no_y_3d(1,:));

c_yx =(k_y.*del_y).*local_dens_no_m_3d(2,:);
c_ym = (k_y*del_y).*local_dens_no_x_3d(2,:)+k_m*k_eff.*Del_col'.*local_dens_no_x_3d(3,:).*(K_y - local_dens_no_x_3d(2,:)); 

c_mx = k_m.*Del_col'.*local_dens_no_y_3d(3,:); 
c_my = k_m.*Del_col'.*local_dens_no_x_3d(3,:);

%c_xm = ((k_x.*del_x)+k_m*k_eff.*Del_col'.*local_dens_no_y_3d(3,:)).*local_dens_no_y_3d(1,:);
%c_ym = ((k_y.*del_y)+k_m*k_eff.*Del_col'.*local_dens_no_x_3d(3,:)).*local_dens_no_x_3d(2,:); % patches with one plant-one frugivore

%all species present
c_xym = k_x*del_x.*local_dens_3d(1,:)+k_m*k_eff.*Del_col'.*local_dens_3d(3,:).*(K_x - local_dens_3d(1,:));
c_yxm = k_y*del_y*local_dens_3d(2,:)+k_m*k_eff.*Del_col'.*local_dens_3d(3,:).*(K_y - local_dens_3d(2,:));
%c_xym = (k_x.*del_x+k_m.*k_eff.*Del_col'.*local_dens_3d(3,:)).*local_dens_3d(1,:); 
%c_yxm = (k_y.*del_y+k_m.*k_eff.*Del_col'.*local_dens_3d(3,:)).*local_dens_3d(2,:);
c_mxy = k_m.*Del_col'.*local_dens_3d(3,:);%all species present


%lambda = c_mx-c_x;
%%
frac_occup_init = [0.1,0.1,0.1];
IC_betweenpatch = repelem(frac_occup_init,num_combinations);
IC_betweenpatch = reshape(IC_betweenpatch,[],num_combinations);
options2 = odeset('NonNegative',[1,2,3]);
if any(e_mxy > 10^4) || any(e_mx > 10^4) || any(e_my > 10^4) || any(e_xym > 10^4) || any(e_yxm > 10^4) 
    c_mx(find((e_mx > 10^4))) = 0; e_mx(find((e_mx > 10^4))) = 0; 
    c_my(find((e_my > 10^4))) = 0; e_my(find((e_my > 10^4))) = 0; 
    c_mxy(find((e_mxy > 10^4))) = 0; e_mxy(find((e_mxy > 10^4))) = 0; 
    c_xym(find((e_xym > 10^4))) = 0; e_xym(find((e_xym > 10^4))) = 0; 
    c_yxm(find((e_yxm > 10^4))) = 0; e_yxm(find((e_yxm > 10^4))) = 0;
end
[t_syst, frac_occup] = ode45(@(t,y)vectorized_BetweenPatchDynamics_allcombos(t,y, c_x0, c_xm, c_xym, c_xy, c_y0, c_ym, c_yxm, c_yx, c_mx, c_mxy, c_my, e_x0, e_xy, e_xm, e_xym, e_y0, e_yx, e_ym, e_yxm, ...
    e_mx, e_mxy, e_my, num_combinations), tspan_meta, IC_betweenpatch(:), options2);

%%
num_tpts_t_syst = length(t_syst);

frac_occup_3d = [];
j = 1;
for i = 1:size(all_params,2)
    %matching_index = find(ismember(all_params(:,i) == Del_and_Q(:,j))
    if good_bad(i) == 1 %3 species coexistence
        %if any(matching_index ~= 0)
        %disp(all_params(:,i))
        %disp('boo-yah!')
        frac_occup_3d = [frac_occup_3d,[frac_occup(end, j);frac_occup(end, j+1);frac_occup(end, j+2)]];
        j = j+3;
    elseif good_bad(i) == -1 % x goes extinct
        frac_occup_3d = [frac_occup_3d,[-1;-1;-1]];
    elseif good_bad(i) == -2
        frac_occup_3d = [frac_occup_3d,[-2;-2;-2]]; %y goes extinct
    elseif good_bad(i) == -3
        frac_occup_3d = [frac_occup_3d,[-3;-3;-3]]; %m goes extinct
    else %more than 1 species going extinct
        %disp('boo-nah?')
        frac_occup_3d = [frac_occup_3d,[-5;-5;-5]];
    end
end
        

%frac_occup_3d = reshape(frac_occup(end,:), [],num_combinations);


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
%save (sprintf('generalist_occupancy_qlow_%d_qhi_%d_delmlo_%d_delmhi_%d_f_%d.mat',q(1)*100,q(end)*100,del_m(1),del_m(end),f*100),'Del','Q','frac_occup_3d','q','del_m','f');
save (sprintf('generalist_occupancy_f_%0.2f.mat',f),'del_m_all','q_all','frac_occup_3d','f','good_bad');


%% surface plot
figure()
surf(del_m_all, q_all, reshape(frac_occup_3d(1,:),numel(q_all),numel(del_m_all)),'EdgeColor','none')
colormap parula; freezeColors;
hold on
surf(del_m_all, q_all, reshape(frac_occup_3d(2,:),numel(q_all),numel(del_m_all)),'EdgeColor','none')
surf(del_m_all, q_all, reshape(frac_occup_3d(3,:),numel(q_all),numel(del_m_all)),'EdgeColor','none')
view(2)
colormap winter; freezeColors;
%surf(Del, Q,occupancy_del_m(:,:,1))
xlabel('mutualist dispersal rate (\delta_m)')
ylabel ('consumption fraction (q)')
%xlim([1.0 29.0]);
zlabel('fraction of patches occupied')
%xline([1.0 2.6, 10.4, 26.9],'--',{'Exploitative (x extinct)','Mutualism (y fitter)','Mutualism (x fitter)', 'Mutualism (y fitter)'})
title(sprintf('Patch occupancy vs mutualist dispersal predation rate; f = %0.2f',f))
legend('Species with mutualist (x)')%, 'Species without mutualist (y)', 'mutualist (m)', 'location', 'best' )
fig1name = sprintf('generalist_vector_occup_q_vs_del_m_f_%0.2f.jpg',f);
print(fig1name,'-djpeg','-r600')
