% These values are to be set. I have set them to arbitrary values

dist_tooltip = 70; %This is only the x distance


curve_point = 50;

gamma_y = 0.2;
gamma_z = 0.2;


r_curve = [curve_point; 0; 0];
curve_length_x = dist_tooltip - r_curve(1);

R_y_c = [cos(gamma_y), 0, sin(gamma_y); 0, 1, 0; -sin(gamma_y), 0, cos(gamma_y)];
R_z_c = [cos(gamma_z), -sin(gamma_z), 0; sin(gamma_z), cos(gamma_z), 0; 0, 0, 1];
R_c = R_z_c * R_y_c;

r_ee_dir = R * [1;0;0];
r_ee_dir_x = r_ee_dir(1);
r_ee_dir = r_ee_dir / r_ee_dir_x * curve_length_x;

r_ee = r_ee_dir + r_curve;

for i=1:size(taskspace,2)
    alpha_y = taskspace(5,i);
    alpha_z = taskspace(6,i);
    
    R_y = [cos(alpha_y), 0, sin(alpha_y); 0, 1, 0; -sin(alpha_y), 0, cos(alpha_y)];
    R_z = [cos(alpha_z), -sin(alpha_z), 0; sin(alpha_z), cos(alpha_z), 0; 0, 0, 1];
    R = R_z * R_y;
    
    temp = R * [1;0;0];
    temp_x = temp(1);
    temp = temp/ temp_x * curve_length_x;
    
    taskspace_t1(1:3,i) = -temp + taskspace(1:3,i);
    
    beta_y = alpha_y - gamma_y;
    beta_z = alpha_z - gamma_z;
    
    R_y_2 = [cos(beta_y), 0, sin(beta_y); 0, 1, 0; -sin(beta_y), 0, cos(beta_y)];
    R_z_2 = [cos(beta_z), -sin(beta_z), 0; sin(beta_z), cos(beta_z), 0; 0, 0, 1];
    R_2 = R_z_2 * R_y_2;
    
    temp2 = R_2 * r_curve;
    taskspace_t2(1:3,i) = -temp2 + taskspace_t1(1:3,i);
    taskspace_t2(4:6,i) = [0; beta_y; beta_z];
   
end
