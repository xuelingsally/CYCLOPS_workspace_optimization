% Euler Angles
ea(1:6) = eaB(1:6);

% Attachment Points
a(:,1) = [eaB(7), radius_tool * sin(ea(1)), radius_tool * cos(ea(1))];
a(:,2) = [eaB(8), radius_tool * sin(ea(2)), radius_tool * cos(ea(2))];
a(:,3) = [eaB(9), radius_tool * sin(ea(3)), radius_tool * cos(ea(3))];
a(:,4) = [eaB(10), radius_tool * sin(ea(4)), radius_tool * cos(ea(4))];
a(:,5) = [eaB(11), radius_tool * sin(ea(5)), radius_tool * cos(ea(5))];
a(:,6) = [eaB(12), radius_tool * sin(ea(6)), radius_tool * cos(ea(6))];

% Base Frame 
B(:,1) = [eaB(13), radius_scaffold * sin(ea(1)), radius_scaffold * cos(ea(1))];
B(:,2) = [eaB(13), radius_scaffold * sin(ea(2)), radius_scaffold * cos(ea(2))];
B(:,3) = [eaB(13), radius_scaffold * sin(ea(3)), radius_scaffold * cos(ea(3))];
B(:,4) = [eaB(14), radius_scaffold * sin(ea(4)), radius_scaffold * cos(ea(4))];
B(:,5) = [eaB(14), radius_scaffold * sin(ea(5)), radius_scaffold * cos(ea(5))];
B(:,6) = [eaB(14), radius_scaffold * sin(ea(6)), radius_scaffold * cos(ea(6))];

gamma_y = eaB(16);
gamma_z = eaB(17);

R_y = [cos(gamma_y), 0, sin(gamma_y); 0, 1, 0; -sin(gamma_y), 0, cos(gamma_y)];
R_z = [cos(gamma_z), -sin(gamma_z), 0; sin(gamma_z), cos(gamma_z), 0; 0, 0, 1];
R = R_z * R_y;

curve_length_x = eaB(15) - eaB(18);

r_ee_dir = R * [1;0;0];
r_ee_dir_x = r_ee_dir(1);
r_ee_dir = r_ee_dir / r_ee_dir_x * curve_length_x;

r_ee = r_ee_dir + [eaB(18); 0; 0];

curve_x = eaB(18);

t_min = ones(6,1) * 1;
t_max = ones(6,1) * 60;

W = zeros(1,6);

f_ee = [0,0,-0.1];
feas = [];
unfeas = [];

for i=1:size(taskspace3,2);
    P = taskspace3(:,i);
    P(1:3,:) = P(1:3,:)/1000;
    P(4:5) = P(5:6);
    P(6) = [];
    
    P = P';
    
    [~, feasible] = feasible_pose(P, a/1000, B/1000, W, f_ee, r_ee'/1000, t_min, t_max);
    
    if feasible == 1
        feas(:,end+1) = P;
    else
        unfeas(:,end+1) = P;
    end
end