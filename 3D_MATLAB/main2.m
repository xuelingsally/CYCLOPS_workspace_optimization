% Main Script
clear;
clc;
%close;

% Constants (All units in mm)
% Radius of overtube/tool
radius_tool = 1;
% Radius of Scaffold
radius_scaffold = 15;
% Length of the scaffold
length_scaffold = 30;
% Lenght of overtube
length_overtube = 15;
% Distance of back of tool to CG of tool
dist_tool_b_cg = 7.5;

% Distance of tool tip to CG of Tool
dist_tooltip = 12;

% Euler Angles
ea(1) = 155/180 * pi;
ea(2) = 280/180 * pi;
ea(3) = 36/180 * pi;
ea(4) = 348/180 * pi;
ea(5) = 242/180 * pi;
ea(6) = 106/180 * pi;

% Attachment Points
a(:,1) = [7.5, radius_tool * cos(ea(1)), radius_tool * sin(ea(1))];
a(:,2) = [7.5, radius_tool * cos(ea(2)), radius_tool * sin(ea(2))];
a(:,3) = [7.5, radius_tool * cos(ea(3)), radius_tool * sin(ea(3))];
a(:,4) = [-7.5, radius_tool * cos(ea(4)), radius_tool * sin(ea(4))];
a(:,5) = [-7.5, radius_tool * cos(ea(5)), radius_tool * sin(ea(5))];
a(:,6) = [-7.5, radius_tool * cos(ea(6)), radius_tool * sin(ea(6))];

% Base Frame 
B(:,1) = [-30, radius_scaffold * cos(ea(1)), radius_scaffold * sin(ea(1))];
B(:,2) = [-30, radius_scaffold * cos(ea(2)), radius_scaffold * sin(ea(2))];
B(:,3) = [-30, radius_scaffold * cos(ea(3)), radius_scaffold * sin(ea(3))];
B(:,4) = [0, radius_scaffold * cos(ea(4)), radius_scaffold * sin(ea(4))];
B(:,5) = [0, radius_scaffold * cos(ea(5)), radius_scaffold * sin(ea(5))];
B(:,6) = [0, radius_scaffold * cos(ea(6)), radius_scaffold * sin(ea(6))];

% Taskspace Definition: List of points that form the boundary of the space
% that the surgeon needs the tool to move in for the operation.
taskspace = [];
% taskspace(:,end+1) = [8, 0, 0];
taskspace(:,end+1) = [5, 1.5, 1.5];
taskspace(:,end+1) = [5, -1.5, 1.5];
taskspace(:,end+1) = [5, -1.5, -1.5];
taskspace(:,end+1) = [5, 1.5, -1.5];
taskspace(:,end+1) = [7.5, 2, -2];
taskspace(:,end+1) = [7.5, 2, 2];
taskspace(:,end+1) = [7.5, -2, -2];
taskspace(:,end+1) = [7.5, -2, 2];
taskspace(:,end+1) = [10, 1.5, -1.5];
taskspace(:,end+1) = [10, 1.5, 1.5];
taskspace(:,end+1) = [10, -1.5, -1.5];
taskspace(:,end+1) = [10, -1.5, 1.5];

draw_cyclops_full(a,B, taskspace, radius_tool, radius_scaffold, length_scaffold, length_overtube, dist_tool_b_cg, dist_tooltip);

%Tensions
t_min = 5 * ones(6,1);
t_max = 60 * ones(6,1);

% Orientation limits
phi_min = [-10/180 * pi; -10/180 * pi];
phi_max = [10/180 * pi; 10/180 * pi];

% Wrench Vector
W = [0,0,0,0,0,0];
r_ee = [dist_tooltip, 0, 0];
%f_ee = [0,0,0];
f_ee = [-1,0,0;
         1,0,0;
        0,-1,0;
         0,1,0;
        0,0,-1;
         0,0,1];

% Workspace Calculation
[wp_size, feasible, unfeasible, t] = dex_workspace(a/1000, B/1000, W, f_ee, r_ee/1000, phi_min, phi_max, t_min, t_max);

% plot
figure;
hold;

plot3(feasible(1,:), feasible(2,:), feasible(3,:), 'bo');
plot3(unfeasible(1,:), unfeasible(2,:), unfeasible(3,:), 'r.');

axis([(min(unfeasible(1,:)) - 20/1000) (max(unfeasible(1,:)) + 20/1000) -(radius_scaffold+5)/1000 (radius_scaffold+5)/1000 -(radius_scaffold+5)/1000 (radius_scaffold+5)/1000]);