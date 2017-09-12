addpath data;
load('./data/L.mat');
load('./data/R.mat');
relative_tp;
radius_tool = 1.75;
radius_scaffold = 26.5491;

figure;
hold;
%% Right tool

taskspace = [];
taskspace = dataR_r';
taskspace(1,:) = taskspace(1,:) + abs(min(taskspace(1,:))) + 5;
%taskspace(2,:) = taskspace(2,:) + radius_scaffold/2 - 7.4;
taskspace(2,:) = taskspace(2,:) + radius_scaffold/2;
taskspace(3,:) = taskspace(3,:) + abs(min(taskspace(3,:))) - radius_scaffold - 3;

eaB = [2.72166; 1.0472; 5.19958; 5.23599; 3.12778; 1.0472; -3.47664; -0.362928; -0.0428612; 13.446; 9.23136; 8.45369; -71.6096; 0; 55.7045; -0.571885; -0.354678; 35.6096];

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

x_middle = (min(B(1,:)) + max(B(1,:)))/2;
p = [x_middle; 0; 0];

%length_scaffold = abs(eaB(13) - eaB(14));
length_scaffold = max(B(1,:)) - min(B(1,:));
length_overtube = max(a(1,:)) - min(a(1,:));

dist_tool_b_cg = abs(min(a(1,:)));



% Draw Overtube;
[y,z,x] = cylinder(radius_tool, 20);
x = x * length_overtube + x_middle + min(a(1,:));
y = y(:,:) - radius_scaffold/2;
s = surf(x,y,z, 'EdgeColor', 'none', 'FaceColor','[0.5,0.5,0.5]');

% Draw tool to curve point
[y,z,x] = cylinder(radius_tool, 20);
x = x * (curve_x - max(a(1,:)));
x = x + max(a(1,:)) + x_middle;
y = y(:,:) - radius_scaffold/2;
s = surf(x,y,z, 'EdgeColor', 'none', 'FaceColor','[0.7,0.7,0.7]');

% Draw tool from curve point to tip
[y,z,x] = cylinder(radius_tool, 20);
circ1(1,:) = x(1,:);
circ1(2,:) = y(1,:);
circ1(3,:) = z(1,:);
circ2(1,:) = x(2,:);
circ2(2,:) = y(2,:);
circ2(3,:) = z(2,:);
for i=1:size(circ1,2)
    circ1(:,i) = R*circ1(:,i) + [curve_x+x_middle;0;0];
    circ2(:,i) = R*circ2(:,i) + r_ee + [x_middle;0;0];
end
x(1,:) = circ1(1,:);
x(2,:) = circ2(1,:);
y(1,:) = circ1(2,:);
y(2,:) = circ2(2,:);
z(1,:) = circ1(3,:);
z(2,:) = circ2(3,:);
y = y(:,:) - radius_scaffold/2;
s = surf(x,y,z, 'EdgeColor', 'none', 'FaceColor','[0.7,0.7,0.7]');



% Draw Tendons
for i=1:size(a,2)
    a_b = a(:,i) + p;
    b = B(:,i);
    
    x = [a_b(1), b(1)];
    y = [a_b(2), b(2)];
    z = [a_b(3), b(3)];
    
    y = y(:,:) - radius_scaffold/2;
    
    plot3(x, y, z, 'r-');
    plot3(x, y, z, 'rx');
end

taskspace(2,:) = taskspace(2,:) - radius_scaffold/2;

% Draw Taskspace
if size(taskspace > 0)
    plot3(taskspace(1,:), taskspace(2,:), taskspace(3,:), '.', 'MarkerEdgeColor', '[0,0.3,0.6]');
end






%% Left tool
taskspace = dataL_r';
taskspace(1,:) = taskspace(1,:) + abs(min(taskspace(1,:))) + 5;
%taskspace(2,:) = taskspace(2,:) + radius_scaffold/2 - 7.4;
taskspace(2,:) = taskspace(2,:) - radius_scaffold/2;
taskspace(3,:) = taskspace(3,:) + abs(min(taskspace(3,:))) - radius_scaffold - 3;

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

x_middle = (min(B(1,:)) + max(B(1,:)))/2;
p = [x_middle; 0; 0];

%length_scaffold = abs(eaB(13) - eaB(14));
length_scaffold = max(B(1,:)) - min(B(1,:));
length_overtube = max(a(1,:)) - min(a(1,:));

dist_tool_b_cg = abs(min(a(1,:)));



% Draw Overtube;
[y,z,x] = cylinder(radius_tool, 20);
x = x * length_overtube + x_middle + min(a(1,:));
y = y(:,:) + radius_scaffold/2;
s = surf(x,y,z, 'EdgeColor', 'none', 'FaceColor','[0.5,0.5,0.5]');

% Draw tool to curve point
[y,z,x] = cylinder(radius_tool, 20);
x = x * (curve_x - max(a(1,:)));
x = x + max(a(1,:)) + x_middle;
y = y(:,:) + radius_scaffold/2;
s = surf(x,y,z, 'EdgeColor', 'none', 'FaceColor','[0.7,0.7,0.7]');

% Draw tool from curve point to tip
[y,z,x] = cylinder(radius_tool, 20);
circ1(1,:) = x(1,:);
circ1(2,:) = y(1,:);
circ1(3,:) = z(1,:);
circ2(1,:) = x(2,:);
circ2(2,:) = y(2,:);
circ2(3,:) = z(2,:);
for i=1:size(circ1,2)
    circ1(:,i) = R*circ1(:,i) + [curve_x+x_middle;0;0];
    circ2(:,i) = R*circ2(:,i) + r_ee + [x_middle;0;0];
end
x(1,:) = circ1(1,:);
x(2,:) = circ2(1,:);
y(1,:) = circ1(2,:);
y(2,:) = circ2(2,:);
z(1,:) = circ1(3,:);
z(2,:) = circ2(3,:);
y = y(:,:) + radius_scaffold/2;
s = surf(x,y,z, 'EdgeColor', 'none', 'FaceColor','[0.7,0.7,0.7]');



% Draw Tendons
for i=1:size(a,2)
    a_b = a(:,i) + p;
    b = B(:,i);
    
    x = [a_b(1), b(1)];
    y = [a_b(2), b(2)];
    z = [a_b(3), b(3)];
    
    y = y(:,:) + radius_scaffold/2;
    
    plot3(x, y, z, 'r-');
    plot3(x, y, z, 'rx');
end

taskspace(2,:) = taskspace(2,:) + radius_scaffold/2;

% Draw Taskspace
if size(taskspace > 0)
    plot3(taskspace(1,:), taskspace(2,:), taskspace(3,:), '.', 'MarkerEdgeColor', '[0.6,0.3,0]');
end


%% Draw combined scaffold

y = [];
z = [];
counter = 1;
for i=-12:12
y(1,counter) = radius_scaffold * cos(i/18*pi);
y(1,counter) = y(1,counter) + radius_scaffold/2;
z(1,counter) = radius_scaffold * sin(i/18*pi);
counter = counter + 1;
end;
for j=7:30
y(1,counter) = radius_scaffold * cos(j/18*pi);
y(1,counter) = y(1,counter) - radius_scaffold/2;
z(1,counter) = radius_scaffold * sin(j/18*pi);
counter = counter + 1;
end;
y(2,:) = y(1,:);
z(2,:) = z(1,:);
x = [];
x(1,1:49) = 0;
x(2,1:49) = -80;
surf(x,y,z, 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'FaceColor', '[1,0.5,0.5]');

