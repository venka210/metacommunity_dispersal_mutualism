clear all; close all;

figure(1); clf
hh = gcf;
set(hh,'PaperUnits','centimeters');
set(hh,'Units','centimeters');
set(gcf,'Position',[0.5 2.0 30 10]) %[starting point (to the right), starting point (from bottom of screen), width, height]

t = tiledlayout(1,3);

%title(t,'Coexistence outcomes with and without dispersal-dependence (DD) tradeoffs')
xlabel(t,'time (local dynamics)','FontSize', 16)
ylabel(t,'local patch densities','FontSize', 16)

%Tile 1
nexttile
load('specialist_occupancy_no_m_q_1.00_del_m_max_10.0.mat')
plot(t_patch,local_dens(:,1),'LineWidth',2); hold on; plot(t_patch, local_dens(:,2),'LineWidth',2); plot(t_patch, local_dens(:,3),'--','LineWidth',2)
title ('a) Patch with x and y only', 'FontSize',12)
clear vars
% xlabel('time (local dynamics)')
% ylabel('local patch density')
%title('Local patch population density vs time')
%legend('Species with mutualist (x)', 'Species without mutualist (y)', 'mutualist (m)', 'location', 'best' )

%Tile 2
nexttile
load('specialist_occupancy_q_1.00_del_m_max_10.0.mat')
plot(t_patch_no_y,local_dens_no_y(:,1),'LineWidth',2); hold on; plot(t_patch_no_y, local_dens_no_y(:,2),'LineWidth',2); plot(t_patch_no_y, local_dens_no_y(:,3),'--','LineWidth',2)
title ('b) Patch with x and m only', 'FontSize',12)
clear vars
% xlabel('time (local dynamics)')
% ylabel('local patch density')
%title('Local patch population density vs time')
%legend('Species with mutualist (x)', 'Species without mutualist (y)', 'mutualist (m)', 'location', 'best' )

%Tile 3
nexttile
load('specialist_occupancy_q_1.00_del_m_max_10.0.mat')
plot(t_patch,local_dens(:,1),'LineWidth',2); hold on; plot(t_patch, local_dens(:,2),'LineWidth',2); plot(t_patch, local_dens(:,3),'--','LineWidth',2)
title ('c) Patch with x, y and m coexisting', 'FontSize',12)
clear vars
% xlabel('time (local dynamics)')
% ylabel('local patch density')
%title('Local patch population density vs time')
%legend('Species with mutualist (x)', 'Species without mutualist (y)', 'mutualist (m)', 'location', 'best' )
%print(sprintf('locdens_vs_time_q_%0.2f_del_m_max_%3.1f.jpg',q,del_m(end)),'-djpeg','-r600')

leg = legend('Species with mutualist (x)', 'Species without mutualist (y)', 'mutualist (m)', 'location', 'best' );
leg.Layout.Tile = 'north';
print -djpeg -r600 Figure_3.jpg