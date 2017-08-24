close all;

radius_tool = 1.75;
%radius_scaffold = 27.3586;
radius_scaffold = 28.64;

load('data_r_tp.mat');
taskspace = data_r_tp;
taskspace(1,:) = taskspace(1,:) + abs(min(taskspace(1,:))) + 5;
taskspace(2,:) = taskspace(2,:) + radius_scaffold/2 - 4;
%taskspace(3,:) = taskspace(3,:) + abs(min(taskspace(3,:))) - radius_scaffold - 3;
taskspace(3,:) = taskspace(3,:) - 13;
data1 = taskspace;

data_t1 = [];
R = [];

r_curve = [eaB(18); 0; 0];

gamma_y = eaB(16);
gamma_z = eaB(17);

curve_length_x = eaB(15) - eaB(18);


R_y_c = [cos(gamma_y), 0, sin(gamma_y); 0, 1, 0; -sin(gamma_y), 0, cos(gamma_y)];
R_z_c = [cos(gamma_z), -sin(gamma_z), 0; sin(gamma_z), cos(gamma_z), 0; 0, 0, 1];
R_c = R_z_c * R_y_c;

r_ee_dir = R_c * [1;0;0];
r_ee_dir_x = r_ee_dir(1);
r_ee_dir = r_ee_dir / r_ee_dir_x * curve_length_x;

% r_ee_dir is the vector from the curving point to the end effector of the tool
curve_length = norm(r_ee_dir);

r_ee = r_ee_dir + r_curve;


for i=1:size(data1,2)
    alpha_y = data1(5,i);
    alpha_z = data1(6,i);
    
    R_y = [cos(alpha_y), 0, sin(alpha_y); 0, 1, 0; -sin(alpha_y), 0, cos(alpha_y)];
    R_z = [cos(alpha_z), -sin(alpha_z), 0; sin(alpha_z), cos(alpha_z), 0; 0, 0, 1];
    R = R_z * R_y;
    
    temp = R * [curve_length;0;0];
    
    data_t1(1:3,i) = -temp + data1(1:3,i);
    
    beta_y = alpha_y - gamma_y;
    beta_z = alpha_z - gamma_z;
    
    R_y_2 = [cos(beta_y), 0, sin(beta_y); 0, 1, 0; -sin(beta_y), 0, cos(beta_y)];
    R_z_2 = [cos(beta_z), -sin(beta_z), 0; sin(beta_z), cos(beta_z), 0; 0, 0, 1];
    R_2 = R_z_2 * R_y_2;
    
    temp2 = R_2 * r_curve;
    data_t2(1:3,i) = -temp2 + data_t1(1:3,i);
    data_t2(4:6,i) = [0;beta_y;beta_z];
end

draw_cyclops_curved(eaB, taskspace, radius_tool, radius_scaffold);
%taskspace2 = data_t1;
%plot3(taskspace2(1,:), taskspace2(2,:), taskspace2(3,:), 'r.');
taskspace3 = data_t2;
plot3(taskspace3(1,:), taskspace3(2,:), taskspace3(3,:), 'g.');
