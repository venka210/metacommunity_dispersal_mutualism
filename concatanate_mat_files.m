clear all
matfiles = dir('occupancy*.mat');

for k = 1:length(matfiles)
  baseFileName = matfiles(k).name;
  matData(k) = load(baseFileName);
end

q_local = unique([matData(1).q,matData(2).q,matData(3).q,matData(4).q]);
del_m_local = unique([matData(1).del_m,matData(2).del_m,matData(3).del_m,matData(4).del_m]);
%frac_occupied_master = matData.frac_occup_3d;
[Del_master, Q_master] = meshgrid(del_m_local,q_local);
arrays = {'array1', 'array2','array3', 'array4';...
    matData(1).frac_occup_3d, matData(2).frac_occup_3d, matData(3).frac_occup_3d,  matData(4).frac_occup_3d};

% Find the array with the maximum size using the max function
%[~, idx_row] = max(cellfun(@(x) size(x,1), arrays(:, 2)));
[~, idx_col] = max(cellfun(@(x) size(x,2), arrays(:, 2)));

% Retrieve information about the array with maximum size
% max_rows_array_name = arrays{idx_row, 1};
% max_rows_array_size = numel(arrays{idx_row, 2});
max_cols_array_name = arrays{idx_col, 1};
max_cols_array_size = size(max_cols_array_name,2);

% %pad arrays such that they all end up to be the same size
% for j = 1:length(matData)
%     if size(matData(j).frac_occup_3d,2) < max_cols_array_size
%         %padarray(matData(j).frac_occup_3d,[3, max_cols_array_size-size(matData(j).frac_occup_3d,2)],0);
%         %matData(j).frac_occup_3d = [matData(j).frac_occup_3d zeros(3,(max_cols_array_size-size(matData(j).frac_occup_3d,2)))];
%     end
% end

%%
figure()
surf(Del_master,Q_master, reshape(frac_occupied_master(1,:),numel(q_local),numel(del_m_local)))
colormap parula; freezeColors;
hold on
surf(Del_master,Q_master, reshape(frac_occupied_master(2,:),numel(q_local),numel(del_m_local)))
colormap winter; freezeColors;
%surf(Del, Q,occupancy_del_m(:,:,1))
xlabel('mutualist dispersal rate (\delta_m)')
ylabel ('consumption fraction (q)')
%xlim([1.0 29.0]);
zlabel('fraction of patches occupied')
%xline([1.0 2.6, 10.4, 26.9],'--',{'Exploitative (x extinct)','Mutualism (y fitter)','Mutualism (x fitter)', 'Mutualism (y fitter)'})
title('Fraction of patches occupied vs mutualist dispersal rate and predation rate')
legend('Species with mutualist (x)')%, 'Species without mutualist (y)', 'mutualist (m)', 'location', 'best' )


    

