load('data_r_tp.mat');
taskspace = data_r_tp;
taskspace(1,:) = taskspace(1,:) + abs(min(taskspace(1,:))) + 5;
taskspace(2,:) = taskspace(2,:) + radius_scaffold/2;
taskspace(3,:) = taskspace(3,:) + abs(min(taskspace(3,:))) - radius_scaffold - 3;
data1 = taskspace;

data_t1 = [];
R = [];

r_curve = [eaB(18); 0; 0];

gamma_y = eaB(15);
gamma_z = eaB(16);

curve_length_x = eaB(15) - eaB(18);

for i=1:size(data1,2)
    phi_y = data1(5,i);
    phi_z = data1(6,i);
    
    R_y = [cos(phi_y), 0, sin(phi_y); 0, 1, 0; -sin(phi_y), 0, cos(phi_y)];
    R_z = [cos(phi_z), -sin(phi_z), 0; sin(phi_z), cos(phi_z), 0; 0, 0, 1];
    R = R_z * R_y;
    
    temp = -R * [1;0;0];
    temp_x = temp(1);
    temp = temp/ temp_x * curve_length_x;
    
    data_t1(1:3,i) = temp + data1(1:3,i);
    
    
    
end

%draw_cyclops_full(a,B, taskspace, radius_tool, radius_scaffold, length_scaffold, length_overtube, dist_tool_b_cg, dist_tooltip);
draw_cyclops_curved(eaB, taskspace, radius_tool, radius_scaffold);
taskspace2 = data_t1;
plot3(taskspace2(1,:), taskspace2(2,:), taskspace2(3,:), 'r.');
