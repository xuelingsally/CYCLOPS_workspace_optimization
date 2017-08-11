% Main Script
clear;
clc;
%close;

% Constants (All units in mm)
% Radius of overtube/tool
radius_tool = 1.75;
% Radius of Scaffold
radius_scaffold = 30;
% Length of the scaffold
length_scaffold = 70;
% Lenght of overtube
length_overtube = 60;
% Distance of back of tool to CG of tool
dist_tool_b_cg = 30;

% Distance of tool tip to CG of Tool
dist_tooltip = 70;

% Euler Angles
ea(1) = 60/180 * pi;
ea(2) = 180/180 * pi;
ea(3) = 300/180 * pi;
ea(4) = 60/180 * pi;
ea(5) = 180/180 * pi;
ea(6) = 300/180 * pi;

% Attachment Points
a(:,1) = [30, radius_tool * cos(ea(1)), radius_tool * sin(ea(1))];
a(:,2) = [30, radius_tool * cos(ea(2)), radius_tool * sin(ea(2))];
a(:,3) = [30, radius_tool * cos(ea(3)), radius_tool * sin(ea(3))];
a(:,4) = [-30, radius_tool * cos(ea(4)), radius_tool * sin(ea(4))];
a(:,5) = [-30, radius_tool * cos(ea(5)), radius_tool * sin(ea(5))];
a(:,6) = [-30, radius_tool * cos(ea(6)), radius_tool * sin(ea(6))];

% Base Frame 
B(:,1) = [0, radius_scaffold * cos(ea(1)), radius_scaffold * sin(ea(1))];
B(:,2) = [0, radius_scaffold * cos(ea(2)), radius_scaffold * sin(ea(2))];
B(:,3) = [0, radius_scaffold * cos(ea(3)), radius_scaffold * sin(ea(3))];
B(:,4) = [-70, radius_scaffold * cos(ea(4)), radius_scaffold * sin(ea(4))];
B(:,5) = [-70, radius_scaffold * cos(ea(5)), radius_scaffold * sin(ea(5))];
B(:,6) = [-70, radius_scaffold * cos(ea(6)), radius_scaffold * sin(ea(6))];

% Taskspace Definition: List of points that form the boundary of the space
% that the surgeon needs the tool to move in for the operation.
taskspace = [];
taskspace(:,end+1) = [20, 5, 5];
taskspace(:,end+1) = [20, -5, 5];
taskspace(:,end+1) = [20, -5, -5];
taskspace(:,end+1) = [20, 5, -5];
taskspace(:,end+1) = [25, 5, -6];
taskspace(:,end+1) = [25, 5, 6];
taskspace(:,end+1) = [25, -5, -6];
taskspace(:,end+1) = [25, -5, 6];
taskspace(:,end+1) = [35, 5, -10];
taskspace(:,end+1) = [30, 5, 5];
taskspace(:,end+1) = [30, -5, -5];
taskspace(:,end+1) = [30, -5, 5];

draw_cyclops_full(a,B, taskspace, radius_tool, radius_scaffold, length_scaffold, length_overtube, dist_tool_b_cg, dist_tooltip);

%Tensions
t_min = 1 * ones(6,1);
t_max = 60 * ones(6,1);

% Orientation limits
phi_min = [-10/180 * pi; -10/180 * pi];
phi_max = [10/180 * pi; 10/180 * pi];

% Wrench Vector
W = [0,0,0,0,0,0];
r_ee = [dist_tooltip, 0, 0];
f_ee = [0,0,0];

% Workspace Calculation
[wp_size, feasible, unfeasible, t] = dex_workspace(a/1000, B/1000, W, f_ee, r_ee/1000, phi_min, phi_max, t_min, t_max, length_scaffold/1000);

plot3(feasible(1,:)*1000, feasible(2,:)*1000, feasible(3,:)*1000, 'b.');

% plot
figure;
hold;

plot3(feasible(1,:), feasible(2,:), feasible(3,:), 'bo');
plot3(unfeasible(1,:), unfeasible(2,:), unfeasible(3,:), 'r.');

axis([-(length_scaffold + 5)/1000 5/1000 -(radius_scaffold+5)/1000 (radius_scaffold+5)/1000 -(radius_scaffold+5)/1000 (radius_scaffold+5)/1000]);