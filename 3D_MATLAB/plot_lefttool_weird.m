eaB = [-1.12297; 2.08043; -2.09433; -2.08292; -0.0122002; 2.04397; -14.0156; -24.0649; -0.0986144; 0.539423; 0; 0.287158; -80; 0; 30.3056; -0.647377; 1.01383; 21.0877];
% Euler Angles
ea(1:6) = eaB(1:6);

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

 
test2a
% Euler Angles
ea(1:6) = eaB(1:6);

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
p = data_t2(1:3,64);
phi_y = data_t2(5,64);
phi_z = data_t2(6,64);
R_y = [cos(phi_y), 0, sin(phi_y); 0, 1, 0; -sin(phi_y), 0, cos(phi_y)];
R_z = [cos(phi_z), -sin(phi_z), 0; sin(phi_z), cos(phi_z), 0; 0, 0, 1];
T_r_overtube = R_z * R_y;


figure;
hold;

length_scaffold = max(B(1,:)) - min(B(1,:));
length_overtube = max(a(1,:)) - min(a(1,:));
[y,z,x] = cylinder(radius_scaffold, 20);
x = x * -length_scaffold + max(B(1,:));
surf(x,y,z, 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'FaceColor', '[1,0.5,0.5]');
[y,z,x] = cylinder(radius_tool, 20);
circ1(1,:) = x(1,:);
circ1(2,:) = y(1,:);
circ1(3,:) = z(1,:);
circ2(1,:) = x(2,:);
circ2(2,:) = y(2,:);
circ2(3,:) = z(2,:);

circ1(1,:) = min(a(1,:)); 
circ2(1,:) = max(a(1,:));
for i=1:size(circ1,2)
circ1(:,i) = T_r_overtube*circ1(:,i) + p;
circ2(:,i) = T_r_overtube*circ2(:,i) + p;
end
x(1,:) = circ1(1,:);
x(2,:) = circ2(1,:);
y(1,:) = circ1(2,:);
y(2,:) = circ2(2,:);
z(1,:) = circ1(3,:);
z(2,:) = circ2(3,:);
s = surf(x,y,z, 'EdgeColor', 'none', 'FaceColor','[0.5,0.5,0.5]');
for i=1:size(a,2)
a_b = T_r_overtube * a(:,i) + p;
b = B(:,i);
    
    x = [a_b(1), b(1)];
    y = [a_b(2), b(2)];
    z = [a_b(3), b(3)];
    
    %y = y(:,:) - radius_scaffold/2;
    
    plot3(x, y, z, 'r-');
    plot3(x, y, z, 'rx');
end;

plot3(taskspace(1,:), taskspace(2,:), taskspace(3,:), 'b.');

x = [taskspace(1,64), data_t1(1,64), data_t2(1,64)];
y = [taskspace(2,64), data_t1(2,64), data_t2(2,64)];
z = [taskspace(3,64), data_t1(3,64), data_t2(3,64)];
plot3(x,y,z, 'k-');

