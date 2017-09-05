function draw_cyclops_curved2(eaB, taskspace, radius_tool, radius_scaffold)

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


add_cables = (size(eaB,1) - 21) / 3;

for i=1:add_cables
    ea(end+1) = eaB(19+(i-1)*3);
    a(:,end+1) = [eaB(20+(i-1)*3), radius_tool * cos(ea(end)), radius_tool * sin(ea(end))];
    B(:,end+1) = [eaB(21+(i-1)*3), radius_scaffold * cos(ea(end)), radius_scaffold * sin(ea(end))];
end

gamma_y = eaB(16);
gamma_z = eaB(17);

R_y = [cos(gamma_y), 0, sin(gamma_y); 0, 1, 0; -sin(gamma_y), 0, cos(gamma_y)];
R_z = [cos(gamma_z), -sin(gamma_z), 0; sin(gamma_z), cos(gamma_z), 0; 0, 0, 1];
R = R_z * R_y;

curve_length_x = eaB(21) - eaB(18);

r_ee_dir = R * [1;0;0];
r_ee_dir_x = r_ee_dir(1);
r_ee_dir = r_ee_dir / r_ee_dir_x * curve_length_x;


gamma2_y = eaB(19);
gamma2_z = eaB(20);

curve_length_x2 = eaB(15) - eaB(21);

R_y_c = [cos(gamma2_y), 0, sin(gamma2_y); 0, 1, 0; -sin(gamma2_y), 0, cos(gamma2_y)];
R_z_c = [cos(gamma2_z), -sin(gamma2_z), 0; sin(gamma2_z), cos(gamma2_z), 0; 0, 0, 1];
R_c = R_z_c * R_y_c;

r_ee_dir2 = R_c * r_ee_dir;

r_ee_dir_x2 = r_ee_dir2(1);
r_ee_dir2 = r_ee_dir2 / r_ee_dir_x2 * curve_length_x2;
% r_ee_dir is the vector from the curving point to the end effector of the tool

r_c2 = r_ee_dir + [eaB(18); 0; 0];

r_ee = r_ee_dir2 + r_ee_dir + [eaB(18); 0; 0];

curve_x = eaB(18);
curve_x2 = eaB(21);

x_middle = (min(B(1,:)) + max(B(1,:)))/2;
p = [x_middle; 0; 0];

length_scaffold = abs(eaB(13) - eaB(14));
length_overtube = max(a(1,:)) - min(a(1,:));

dist_tool_b_cg = abs(min(a(1,:)));




figure;
hold;

% Draw scaffold;
[y,z,x] = cylinder(radius_scaffold, 20);
x = x * -length_scaffold + max(B(1,:));
surf(x,y,z, 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'FaceColor', '[1,0.5,0.5]');

% Draw Overtube;
[y,z,x] = cylinder(radius_tool, 20);
x = x * length_overtube + x_middle + min(a(1,:));
s = surf(x,y,z, 'EdgeColor', 'none', 'FaceColor','[0.5,0.5,0.5]');

% Draw tool to curve point
% [y,z,x] = cylinder(radius_tool, 20);
% x = x * (curve_x - max(a(1,:)));
% x = x + max(a(1,:)) + x_middle;
% s = surf(x,y,z, 'EdgeColor', 'none', 'FaceColor','[0.7,0.7,0.7]');

x = [x_middle,curve_x + x_middle];
y = [0, 0];
z = [0,0];
plot3(x, y, z, 'k-', 'LineWidth', 1);

% % Draw tool from curve point to tip
% [y,z,x] = cylinder(radius_tool, 20);
% circ1(1,:) = x(1,:);
% circ1(2,:) = y(1,:);
% circ1(3,:) = z(1,:);
% circ2(1,:) = x(2,:);
% circ2(2,:) = y(2,:);
% circ2(3,:) = z(2,:);
% for i=1:size(circ1,2)
%     circ1(:,i) = R*circ1(:,i) + [curve_x+x_middle;0;0];
%     circ2(:,i) = R*circ2(:,i) + r_ee + [x_middle;0;0];
% end
% x(1,:) = circ1(1,:);
% x(2,:) = circ2(1,:);
% y(1,:) = circ1(2,:);
% y(2,:) = circ2(2,:);
% z(1,:) = circ1(3,:);
% z(2,:) = circ2(3,:);
% s = surf(x,y,z, 'EdgeColor', 'none', 'FaceColor','[0.7,0.7,0.7]');

x = [curve_x + x_middle, r_c2(1)+x_middle];
y = [0, r_c2(2)];
z = [0, r_c2(3)];
plot3(x, y, z, 'k-', 'LineWidth', 1);

x = [r_c2(1)+x_middle, r_ee(1)+x_middle];
y = [r_c2(2), r_ee(2)];
z = [r_c2(3), r_ee(3)];
plot3(x, y, z, 'k-', 'LineWidth', 1);



% Draw Tendons
for i=1:size(a,2)
    a_b = a(:,i) + p;
    b = B(:,i);
    
    x = [a_b(1), b(1)];
    y = [a_b(2), b(2)];
    z = [a_b(3), b(3)];
    
    plot3(x, y, z, 'r-');
    plot3(x, y, z, 'rx');
end

% Draw Taskspace
if size(taskspace > 0)
    plot3(taskspace(1,:), taskspace(2,:), taskspace(3,:), 'b.');
end

max_tp_x = max(taskspace(1,:));
max_x_axis = max([max_tp_x, r_ee(1)+x_middle]);

axis([min(B(1,:))-5 (max_x_axis+10) -(radius_scaffold+5) radius_scaffold+5 -(radius_scaffold+5) radius_scaffold+5 ]);
%axis([-85 50 -(radius_scaffold+5) radius_scaffold+5 -(radius_scaffold+5) radius_scaffold+5 ]);

end