clear all; close all;

figure(1); clf
hh = gcf;
set(hh,'PaperUnits','centimeters');
set(hh,'Units','centimeters');
set(gcf,'Position',[0.5 1 20 15])

t = tiledlayout('flow');

%title(t,'Coexistence outcomes with and without dispersal-dependence (DD) tradeoffs')
xlabel(t,'Dispersal rates of the mutualist (\delta_m)','FontSize', 10)
ylabel(t,'Patch occupancies of species (p_x, p_y, p_m)','FontSize', 10)

%Tile 1
nexttile
load('specialist_occupancy_q_1.00_del_m_max_10.0.mat')
plot(del_m, occupancy_del_m(:, 1:3),'LineWidth',3)
yline(0.1368,'g--','LineWidth',3)
%axis tight
%set(gca,'YDir','normal','FontSize',16)
% xlabel('F_2 competitive ability (\tau_{12})','FontSize', 16);
% ylabel('F_1 competitive ability (\tau_{21})','FontSize', 16);
title ('a) Frugivory shifts from antagonism to mutualism when x is poor disperser (\delta_x = 0.01)','FontSize',8)


legend('Species with mutualist (x)', 'Species without mutualist (y)', 'mutualist (m)', 'p_x w/o mutualist','Location','southeastoutside', 'FontSize', 10)


%Tile 2
nexttile
load('specialist_occupancy_ext_q_1.00_del_m_max_10.0.mat')
plot(del_m, occupancy_del_m(:, 1:3),'LineWidth',3)
yline(0,'g--','LineWidth',3)
%axis tight
%set(gca,'YDir','normal','FontSize',16)
% xlabel('F_2 competitive ability (\tau_{12})','FontSize', 16);
% ylabel('F_1 competitive ability (\tau_{21})','FontSize', 16);
title ('b) Frugivore engenders coexistence when x is a very poor disperser (\delta_x = 0.001)','FontSize',8)

print -djpeg -r600 fig_3_bothcases.jpg