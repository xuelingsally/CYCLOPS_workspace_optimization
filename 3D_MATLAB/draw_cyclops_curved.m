function draw_cyclops_curved(eaB, taskspace, radius_tool, radius_scaffold)

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

r_ee = [eaB(15); eaB(16); eaB(17)];

curve_x = eaB(18);

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
surf(x,y,z, 'FaceAlpha', 0.4, 'EdgeColor', 'none', 'FaceColor', 'y');

% Draw Overtube;
[y,z,x] = cylinder(radius_tool, 20);
x = x * length_overtube + x_middle + min(a(1,:));
s = surf(x,y,z, 'EdgeColor', 'none', 'FaceColor','k');

% Draw tool to curve point
[y,z,x] = cylinder(radius_tool, 20);
temp = (length_overtube + x_middle + min(a(1,:)));
x = temp + x * (curve_x + x_middle - temp);
s = surf(x,y,z, 'EdgeColor', 'none', 'FaceColor','r');

% Draw tool from curve point to tip
x = [curve_x + x_middle, r_ee(1)];
y = [0, r_ee(2)];
z = [0, r_ee(3)];
plot3(x, y, z, 'r-', 'LineWidth', 1);

% Draw Tendons
for i=1:size(a,2)
    a_b = a(:,i) + p;
    b = B(:,i);
    
    x = [a_b(1), b(1)];
    y = [a_b(2), b(2)];
    z = [a_b(3), b(3)];
    
    plot3(x, y, z, 'g-');
end

% Draw Taskspace
plot3(taskspace(1,:), taskspace(2,:), taskspace(3,:), 'b.');

max_tp_x = max(taskspace(1,:));
max_x_axis = max([max_tp_x, r_ee(1)]);

axis([min(B(1,:))-5 (max_x_axis+30) -(radius_scaffold+5) radius_scaffold+5 -(radius_scaffold+5) radius_scaffold+5 ]);

end