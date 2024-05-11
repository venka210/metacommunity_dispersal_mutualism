%%This code produces Fig 4 in main text. If provided .mat files are being used (which has condensed the data), then
%simply run the code.
clear all; close all;

figure(1); clf
hh = gcf;
set(hh,'PaperUnits','centimeters');
set(hh,'Units','centimeters');
set(gcf,'Position',[0.5 3 33 12])

t = tiledlayout(1,2);

%title(t,'Coexistence outcomes with and without dispersal-dependence (DD) tradeoffs')

%Tile 1
nexttile
load('generalist_occupancy_f_1.00.mat')
x_occ = reshape(frac_occup_3d(1,:),numel(q_all),numel(del_m_all));
y_occ = reshape(frac_occup_3d(2,:),numel(q_all),numel(del_m_all));
del_4 = x_occ(:,41)-y_occ(:,41);
del_7 = x_occ(:,71)-y_occ(:,71);
del_10 = x_occ(:,101)-y_occ(:,101);
plot(q_all, del_4, 'LineWidth',2); hold on; plot(q_all, del_7,'LineWidth',2); plot(q_all,del_10, 'LineWidth',2)
yline(0)
xlabel('Consumption fraction (q)')
ylabel('Difference in patch occupancies (p_x - p_y)')
legend('Low \delta_m', 'Intermediate \delta_m', 'High \delta_m','location','best')
title ('a)', 'FontSize',12, 'HorizontalAlignment','left')

%Tile 2
nexttile
surf(del_m_all, q_all, reshape(frac_occup_3d(1,:),numel(q_all),numel(del_m_all)),'EdgeColor','none')
colormap parula; freezeColors;
hold on
surf(del_m_all, q_all, reshape(frac_occup_3d(2,:),numel(q_all),numel(del_m_all)),'EdgeColor','none')
surf(del_m_all, q_all, reshape(frac_occup_3d(3,:),numel(q_all),numel(del_m_all)),'EdgeColor','none')
view(2)
colormap winter; freezeColors;
xlabel('mutualist dispersal rate (\delta_m)')
ylabel ('consumption fraction (q)')
title ('b)', 'FontSize',12,'HorizontalAlignment','left')
%title (title(sprintf('a) Patch occupancy vs mutualist dispersal predation rate; f = %0.2f',f)), 'FontSize',12)

print -djpeg  -r600 Figure_4.jpg