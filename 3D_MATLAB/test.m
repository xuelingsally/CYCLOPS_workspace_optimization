load('relative_data.mat');

data1 = relative_data;
data1(1,:) = data1(1,:) + 18;

data1_t = [];

for i=1:size(data1,2)
%     r_ee_g_z = sin(data1(6,i)) * dist_tooltip;
%     h_y = cos(data1(6,i)) * dist_tooltip;
%     r_ee_g_y = sin(data1(5,i)) * h_y;
%     r_ee_g_x = cos(data1(5,i)) * h_y;
%     data1_t(:,i) = [data1(1,i) - r_ee_g_x; data1(2,i) - r_ee_g_y; data1(3,i) + r_ee_g_z; data1(4,i); data1(5,i); data1(6,i)];
    r_ee_g_y = sin(data1(5,i));
    r_ee_g_z = sin(data1(6,i)) * cos(data1(5,i));
    r_ee_g_x = cos(data1(6,i)) * cos(data1(5,i));
    
    r_ee = [r_ee_g_x; r_ee_g_y; r_ee_g_z];
    r_ee = -r_ee/norm(r_ee) * dist_tooltip;
    
    data1_t(:,i) = [data1(1,i) + r_ee(1); data1(2,i) + r_ee(2); data1(3,i) + r_ee(3); data1(4,i); data1(5,i); data1(6,i)];
end

taskspace = data1;
draw_cyclops_full(a,B, taskspace, radius_tool, radius_scaffold, length_scaffold, length_overtube, dist_tool_b_cg, dist_tooltip);
taskspace2 = data1_t;
plot3(taskspace2(1,:), taskspace2(2,:), taskspace2(3,:), 'r.');
