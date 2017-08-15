load('relative_data.mat');

data1 = relative_data;
data1(1,:) = data1(1,:) + 18;

data1_t = [];

for i=1:size(data1,2)
r_ee_g_z = sin(data1(6,i)) * dist_tooltip;
h_y = cos(data1(6,i)) * dist_tooltip;
r_ee_g_y = sin(data1(5,i)) * h_y;
r_ee_g_x = cos(data1(5,i)) * h_y;
data1_t(:,i) = [data1(1,i) - r_ee_g_x; data1(2,i) - r_ee_g_y; data1(3,i) - r_ee_g_z; data1(4,i); data1(5,i); data1(6,i)];
end;