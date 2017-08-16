load('data_r_tp.mat');
taskspace = data_r_tp;
taskspace(1,:) = taskspace(1,:) + abs(min(taskspace(1,:))) + 5;
% data1 = relative_data;
% data1(1,:) = data1(1,:) + 18;
data1 = taskspace;

data1_t = [];
R = [];

for i=1:size(data1,2)
    phi_y = data1(5,i);
    phi_z = data1(6,i);
    
    R_y = [cos(phi_y), 0, sin(phi_y); 0, 1, 0; -sin(phi_y), 0, cos(phi_y)];
    R_z = [cos(phi_z), -sin(phi_z), 0; sin(phi_z), cos(phi_z), 0; 0, 0, 1];
    R = R_y * R_z;
    
    r_ee_temp = [dist_tooltip;0;0];
    
    data1_t(1:3,i) = -R * r_ee_temp + data1(1:3,i);
    data1_t(4:6,i) = data1(4:6,i);
end

taskspace = data1;
draw_cyclops_full(a,B, taskspace, radius_tool, radius_scaffold, length_scaffold, length_overtube, dist_tool_b_cg, dist_tooltip);
taskspace2 = data1_t;
plot3(taskspace2(1,:), taskspace2(2,:), taskspace2(3,:), 'r.');
