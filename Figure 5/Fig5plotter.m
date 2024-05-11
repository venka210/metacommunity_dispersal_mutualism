%%This code produces Fig 5 in main text. If provided .mat files are being used (which has condensed the data), then
%simply run the code.
clear all; close all;

figure(1); clf
hh = gcf;
set(hh,'PaperUnits','centimeters');
set(hh,'Units','centimeters');
set(gcf,'Position',[0.5 1 43 25])

t = tiledlayout(3,2);

%title(t,'Coexistence outcomes with and without dispersal-dependence (DD) tradeoffs')
xlabel(t,'mutualist dispersal rate (\delta_m)','FontSize', 16)
ylabel(t,'mutualist consumption fraction (q)','FontSize', 16)

%Tile 1
nexttile
load('generalist_occupancy_f_0.08.mat')
surf(del_m_all, q_all, reshape(frac_occup_3d(1,:),numel(q_all),numel(del_m_all)),'EdgeColor','none')
colormap parula; freezeColors;
hold on
surf(del_m_all, q_all, reshape(frac_occup_3d(2,:),numel(q_all),numel(del_m_all)),'EdgeColor','none')
surf(del_m_all, q_all, reshape(frac_occup_3d(3,:),numel(q_all),numel(del_m_all)),'EdgeColor','none')
view(2)
colormap winter; freezeColors;
title ('a) f = 0.08', 'FontSize',12)
clear vars

%Tile 2
nexttile
load('generalist_occupancy_f_0.16.mat')
surf(del_m_all, q_all, reshape(frac_occup_3d(1,:),numel(q_all),numel(del_m_all)),'EdgeColor','none')
colormap parula; freezeColors;
hold on
surf(del_m_all, q_all, reshape(frac_occup_3d(2,:),numel(q_all),numel(del_m_all)),'EdgeColor','none')
surf(del_m_all, q_all, reshape(frac_occup_3d(3,:),numel(q_all),numel(del_m_all)),'EdgeColor','none')
view(2)
colormap winter; freezeColors;
title ('b) f = 0.16', 'FontSize',12)
clear vars

%Tile 3
nexttile
load('generalist_occupancy_f_0.44.mat')
surf(del_m_all, q_all, reshape(frac_occup_3d(1,:),numel(q_all),numel(del_m_all)),'EdgeColor','none')
colormap parula; freezeColors;
hold on
surf(del_m_all, q_all, reshape(frac_occup_3d(2,:),numel(q_all),numel(del_m_all)),'EdgeColor','none')
surf(del_m_all, q_all, reshape(frac_occup_3d(3,:),numel(q_all),numel(del_m_all)),'EdgeColor','none')
view(2)
colormap winter; freezeColors;
title ('c) f = 0.44', 'FontSize',12)
clear vars

%Tile 4
nexttile
load('generalist_occupancy_f_0.48.mat')
surf(del_m_all, q_all, reshape(frac_occup_3d(1,:),numel(q_all),numel(del_m_all)),'EdgeColor','none')
colormap parula; freezeColors;
hold on
surf(del_m_all, q_all, reshape(frac_occup_3d(2,:),numel(q_all),numel(del_m_all)),'EdgeColor','none')
surf(del_m_all, q_all, reshape(frac_occup_3d(3,:),numel(q_all),numel(del_m_all)),'EdgeColor','none')
view(2)
colormap winter; freezeColors;
title ('d) f = 0.48', 'FontSize',12)
clear vars

%Tile 5
nexttile
load('generalist_occupancy_f_0.88.mat')
surf(del_m_all, q_all, reshape(frac_occup_3d(1,:),numel(q_all),numel(del_m_all)),'EdgeColor','none')
colormap parula; freezeColors;
hold on
surf(del_m_all, q_all, reshape(frac_occup_3d(2,:),numel(q_all),numel(del_m_all)),'EdgeColor','none')
surf(del_m_all, q_all, reshape(frac_occup_3d(3,:),numel(q_all),numel(del_m_all)),'EdgeColor','none')
view(2)
colormap winter; freezeColors;
title ('e) f = 0.88', 'FontSize',12)
clear vars

%Tile 6
nexttile
load('generalist_occupancy_f_0.96.mat')
surf(del_m_all, q_all, reshape(frac_occup_3d(1,:),numel(q_all),numel(del_m_all)),'EdgeColor','none')
colormap parula; freezeColors;
hold on
surf(del_m_all, q_all, reshape(frac_occup_3d(2,:),numel(q_all),numel(del_m_all)),'EdgeColor','none')
surf(del_m_all, q_all, reshape(frac_occup_3d(3,:),numel(q_all),numel(del_m_all)),'EdgeColor','none')
view(2)
colormap winter; freezeColors;
title ('f) f = 0.96', 'FontSize',12)
clear vars

print -djpeg -r600 fig_5_generalist.jpg
