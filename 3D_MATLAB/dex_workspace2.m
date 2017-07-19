function [wp_size, feasible, unfeasible, t] = dex_workspace2(a, B, W, f_ee, r_ee, phi_min, phi_max, t_min, t_max)

x_res = 5;
theta_res = 12;
r_res = 5;

% Determine Search Volume
x_middle = (min(B(1,:)) + max(B(1,:)))/2;
p_arbitrary = [x_middle; 0; 0];

x_space_length1 = x_space_length(a, B, p_arbitrary);
x_space_length2 = x_space_length(-a, -B, -p_arbitrary);

x_space = [x_middle - x_space_length1, x_middle + x_space_length2];

x_step = (x_space_length1 + x_space_length2)/x_res;
theta_step = 2 * pi / theta_res;
radius = norm(B(2:3));
r_step = radius/r_res;

x_temp = x_space(1);
theta_temp = 0;
r_temp = 0;

vol_grid = [];
for i=1:x_res + 1
    for j=1:theta_res
        for k=1:r_res + 1
           vol_grid(:,end+1) = [x_temp; r_temp*cos(theta_temp); r_temp*sin(theta_temp)];
           r_temp = r_temp + r_step;
        end
        r_temp = 0;
        theta_temp = theta_temp + theta_step;
    end
    theta_temp = 0;
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


end