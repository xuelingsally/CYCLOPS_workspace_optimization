function val = my_objective_function3(eaB, W, f_ee, r_ee, phi_min, phi_max, t_min, t_max, taskspace, radius_tool, radius_scaffold)

val = 0;

% Euler Angles
ea(1) = eaB(1);
ea(2) = eaB(2);
ea(3) = eaB(3);
ea(4) = eaB(4);
ea(5) = eaB(5);
ea(6) = eaB(6);

% Attachment Points
a(:,1) = [eaB(7), radius_tool * cos(ea(1)), radius_tool * sin(ea(1))];
a(:,2) = [eaB(8), radius_tool * cos(ea(2)), radius_tool * sin(ea(2))];
a(:,3) = [eaB(9), radius_tool * cos(ea(3)), radius_tool * sin(ea(3))];
a(:,4) = [eaB(10), radius_tool * cos(ea(4)), radius_tool * sin(ea(4))];
a(:,5) = [eaB(11), radius_tool * cos(ea(5)), radius_tool * sin(ea(5))];
a(:,6) = [eaB(12), radius_tool * cos(ea(6)), radius_tool * sin(ea(6))];

% Base Frame 
B(:,1) = [eaB(13), radius_scaffold * cos(ea(1)), radius_scaffold * sin(ea(1))];
B(:,2) = [eaB(13), radius_scaffold * cos(ea(2)), radius_scaffold * sin(ea(2))];
B(:,3) = [eaB(13), radius_scaffold * cos(ea(3)), radius_scaffold * sin(ea(3))];
B(:,4) = [eaB(14), radius_scaffold * cos(ea(4)), radius_scaffold * sin(ea(4))];
B(:,5) = [eaB(14), radius_scaffold * cos(ea(5)), radius_scaffold * sin(ea(5))];
B(:,6) = [eaB(14), radius_scaffold * cos(ea(6)), radius_scaffold * sin(ea(6))];


% Check for Crossing of Cables: i.e. a_x_i < a_x_j, B_x_i > B_x_j, B_y_i = B_y_j;
for i=1:size(a,2)
    for j=1:size(a,2)
        if (a(1,i) < a(1,j) && B(1,i) > B(1,j))
            %val = -1000;
            val = -100000 * abs(a(1,j) - a(1,i)) - 100;
            return;
        end
    end
end

%Check if the taskspace is reachable by the tool
feasible_temp = zeros(5,1);
for j=1:size(f_ee, 1)
    for i=1:size(taskspace, 2)
        P = [taskspace(1,i), taskspace(2,i), taskspace(3,i), 0, 0];
        [~, feasible_temp(1)] =  feasible_pose(P, a/1000, B/1000, W, f_ee(j,:), r_ee/1000, t_min, t_max);
        P = [taskspace(1,i), taskspace(2,i), taskspace(3,i), phi_min(1), 0];
        [~, feasible_temp(2)] =  feasible_pose(P, a/1000, B/1000, W, f_ee(j,:), r_ee/1000, t_min, t_max);
        P = [taskspace(1,i), taskspace(2,i), taskspace(3,i), phi_max(1), 0];
        [~, feasible_temp(3)] =  feasible_pose(P, a/1000, B/1000, W, f_ee(j,:), r_ee/1000, t_min, t_max);
        P = [taskspace(1,i), taskspace(2,i), taskspace(3,i), phi_min(2), 0];
        [~, feasible_temp(4)] =  feasible_pose(P, a/1000, B/1000, W, f_ee(j,:), r_ee/1000, t_min, t_max);
        P = [taskspace(1,i), taskspace(2,i), taskspace(3,i), phi_max(2), 0];
        [~, feasible_temp(5)] =  feasible_pose(P, a/1000, B/1000, W, f_ee(j,:), r_ee/1000, t_min, t_max);

        if sum(feasible_temp) ~= 5
            val = val -1;
        end
    end
end

if val < 0
    return
end

[wp_size, ~, ~, ~] = dex_workspace(a/1000, B/1000, W, [0,0,0], r_ee/1000, phi_min, phi_max, t_min, t_max);
val = wp_size;

end