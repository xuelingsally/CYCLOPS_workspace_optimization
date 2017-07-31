function [wp_size, feasible, unfeasible, t] = dex_workspace(a, B, W, f_ee, r_ee, phi_min, phi_max, t_min, t_max, length_scaffold)

x_res = 10;
y_res = 10;
z_res = 10;

% Determine Search Volume
x_middle = (min(B(1,:)) + max(B(1,:)))/2;
p_arbitrary = [x_middle; 0; 0];

x_space_length1 = x_space_length(a, B, p_arbitrary);
x_space_length2 = x_space_length(-a, -B, -p_arbitrary);

x_space = [x_middle - x_space_length1, x_middle + x_space_length2];

radius = norm(B(2:3,1));

x_step = length_scaffold/x_res;
y_step = 2 * radius /y_res;
z_step = 2 * radius /z_res;


x_temp = x_space(1);
y_temp = -radius;
z_temp = -radius;

vol_grid = [];
%for i=1:x_res + 1
while x_temp <= x_space(2); 
    for j=1:y_res + 1
        for k=1:z_res + 1
           if (z_temp * z_temp + y_temp * y_temp) <= (radius * radius)
                vol_grid(:,end+1) = [x_temp; y_temp; z_temp];
           end
           z_temp = z_temp + z_step;
        end
        z_temp = -radius;
        y_temp = y_temp + y_step;
    end
    y_temp = -radius;
    x_temp = x_temp + x_step;
end

feasible_temp = zeros(5 * size(f_ee,1),1);
feasible = [];
unfeasible = [];

for i=1:size(vol_grid,2)
    for k=1:size(f_ee,1)
        P = [vol_grid(1,i), vol_grid(2,i), vol_grid(3,i), 0, 0];
        [t(:,i), feasible_temp(1 + (k - 1) * 5)] = feasible_pose(P, a, B, W, f_ee(k,:), r_ee, t_min, t_max);
        P = [vol_grid(1,i), vol_grid(2,i), vol_grid(3,i), phi_min(1), 0];
        [~, feasible_temp(2 + (k - 1) * 5)] = feasible_pose(P, a, B, W, f_ee(k,:), r_ee, t_min, t_max);
        P = [vol_grid(1,i), vol_grid(2,i), vol_grid(3,i), phi_max(1), 0];
        [~, feasible_temp(3 + (k - 1) * 5)] = feasible_pose(P, a, B, W, f_ee(k,:), r_ee, t_min, t_max);
        P = [vol_grid(1,i), vol_grid(2,i), vol_grid(3,i), 0, phi_min(2)];
        [~, feasible_temp(4 + (k - 1) * 5)] = feasible_pose(P, a, B, W, f_ee(k,:), r_ee, t_min, t_max);
        P = [vol_grid(1,i), vol_grid(2,i), vol_grid(3,i), 0, phi_max(2)];
        [~, feasible_temp(5 + (k - 1) * 5)] = feasible_pose(P, a, B, W, f_ee(k,:), r_ee, t_min, t_max);
    end
    
    if sum(feasible_temp) == 5 * size(f_ee,1)
        feasible(:,end+1) = vol_grid(:,i);
    else
        unfeasible(:,end+1) = vol_grid(:,i);
    end
end

x_search = x_space_length1 + x_space_length2;
yz_search = pi * radius * radius;
search_vol = x_search * yz_search;
wp_size = size(feasible, 2) / size(vol_grid, 2) * search_vol;
% wp_size = size(feasible, 2)/ (2000);

end