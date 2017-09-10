close all;

addpath data;
load('./data/L.mat');
load('./data/R.mat');
relative_tp;

radius_tool = 1.75;
radius_scaffold = 26.5491;
%radius_scaffold = 22.7564;
%radius_scaffold = 28.64;
%radius_scaffold = 30.51;

%load('data_r_tp.mat');
%taskspace = data_r_tp;
taskspace = dataR_r';
taskspace(1,:) = taskspace(1,:) + abs(min(taskspace(1,:))) + 5;
%taskspace(2,:) = taskspace(2,:) + radius_scaffold/2 - 7.4;
taskspace(2,:) = taskspace(2,:) + radius_scaffold/2;
taskspace(3,:) = taskspace(3,:) + abs(min(taskspace(3,:))) - radius_scaffold - 3;
%taskspace(3,:) = taskspace(3,:) - 13;
data1 = taskspace;

data_t1 = [];
R = [];


%Find r_ee
gamma_y = eaB(16);
gamma_z = eaB(17);
r_ee_x = eaB(15);

R_y = [cos(gamma_y), 0, sin(gamma_y); 0, 1, 0; -sin(gamma_y), 0, cos(gamma_y)];
R_z = [cos(gamma_z), -sin(gamma_z), 0; sin(gamma_z), cos(gamma_z), 0; 0, 0, 1];
R = R_z * R_y;

r_ee_dir = R * [1;0;0];
r_ee_dir_x = r_ee_dir(1);
r_ee = r_ee_dir / r_ee_dir_x * r_ee_x;

r_ee_length = norm(r_ee);

beta_y = eaB(18) - gamma_y;
beta_z = eaB(19) - gamma_z;

for i=1:size(data1,2)
    alpha_y = data1(5,i);
    alpha_z = data1(6,i);
    
    
    R_y = [cos(alpha_y - beta_y), 0, sin(alpha_y - beta_y); 0, 1, 0; -sin(alpha_y - beta_y), 0, cos(alpha_y - beta_y)];
    R_z = [cos(alpha_z - beta_z), -sin(alpha_z - beta_z), 0; sin(alpha_z - beta_z), cos(alpha_z - beta_z), 0; 0, 0, 1];
    R = R_z * R_y;
    
    temp = R * [r_ee_length;0;0];
    
    data_t1(1:3,i) = -temp + data1(1:3,i);
    
    data_t1(4:6,i) = [0;alpha_y - eaB(18); alpha_z - eaB(19)];

end


draw_cyclops_curvedG(eaB, taskspace, radius_tool, radius_scaffold);
taskspace2 = data_t1;
%plot3(taskspace2(1,:), taskspace2(2,:), taskspace2(3,:), 'g.');