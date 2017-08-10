function draw_cyclops_full(a,B, taskspace, radius_tool, radius_scaffold, length_scaffold, length_overtube, dist_tool_b_cg, dist_tooltip)

x_middle = (min(B(1,:)) + max(B(1,:)))/2;
p = [x_middle; 0; 0];

figure;
hold;

% Draw scaffold;
[y,z,x] = cylinder(radius_scaffold, 20);
x = x * -length_scaffold;
surf(x,y,z, 'FaceAlpha', 0.4, 'EdgeColor', 'none', 'FaceColor', 'y');

% Draw Overtube;
[y,z,x] = cylinder(radius_tool, 20);
x = x * length_overtube + x_middle - dist_tool_b_cg;
s = surf(x,y,z, 'EdgeColor', 'none', 'FaceColor','k');

% Draw tool;
[y,z,x] = cylinder(radius_tool, 20);
x = x * (dist_tooltip + x_middle - (length_overtube + x_middle - dist_tool_b_cg)) + (length_overtube + x_middle - dist_tool_b_cg);
surf(x,y,z, 'EdgeColor', 'none', 'FaceColor', 'r');


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

%axis([-length_scaffold-5 (dist_tooltip + x_middle + (length_overtube + x_middle - dist_tool_b_cg) +10) -(radius_scaffold+5) radius_scaffold+5 -(radius_scaffold+5) radius_scaffold+5 ]);
axis([-length_scaffold-5 (max(taskspace(1,:))+10) -(radius_scaffold+5) radius_scaffold+5 -(radius_scaffold+5) radius_scaffold+5 ]);

end
