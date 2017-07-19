% Main Script
clear;
clc;
%close;

addpath PSwarmM_v2_1;

%% Constants (All units in mm)
% Radius of overtube/tool
radius_tool = 1;
% Radius of Scaffold
radius_scaffold = 30;
% Length of the scaffold
length_scaffold = 100;
% Lenght of overtube
length_overtube = 60;
% Distance of back of tool to CG of tool
dist_tool_b_cg = 30;

% Distance of tool tip to CG of Tool
% dist_tooltip = 50;


%Tensions
t_min = 5 * ones(6,1);
t_max = 60 * ones(6,1);

% Orientation limits
phi_min = [-10/180 * pi; -10/180 * pi];
phi_max = [10/180 * pi; 10/180 * pi];

% Wrench Vector
W = [0,0,0,0,0,0];
r_ee = [0, 0, 0];
f_ee = [-1,0,0;
         1,0,0;
        0,-1,0;
         0,1,0;
        0,0,-1;
         0,0,1];

%% Taskspace Definition: List of points that form the boundary of the space
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
taskspace(:,end+1) = [30, 5, -5];
taskspace(:,end+1) = [35, 5, -10];
taskspace(:,end+1) = [30, 5, 5];
taskspace(:,end+1) = [30, -5, -5];
taskspace(:,end+1) = [30, -5, 5];

% Based on the tooltip, find the poses that the CG of the overtube+tool is
% required to reach.
% for i=1:size(taskspace,2)
%     taskspace_translated(:,i) = taskspace(:,i) - [dist_tooltip; 0; 0];
% end
% taskspace_d = taskspace_translated;
taskspace_d = taskspace;

%% Design Vector
eaB_min = ones(15,1);
eaB_max = ones(15,1);

%Limits for e
eaB_min(1:6) = eaB_min(1:6) * 0;
eaB_max(1:6) = eaB_max(1:6) * 2 * pi;
%Limits for a
eaB_min(7:9) = eaB_min(7:9) * (-dist_tool_b_cg);
eaB_max(7:9) = eaB_max(7:9) * 0;
% Back 3 Tendons
eaB_min(10:12) = eaB_min(10:12) * 0;
eaB_max(10:12) = eaB_max(10:12) * (length_overtube - dist_tool_b_cg);
%Limits for B
eaB_min(13:14) = eaB_min(13:14) * (-length_scaffold);
eaB_max(13:14) = eaB_max(13:14) * 0;

%Limits for distance tooltip
eaB_min(15) = eaB_min(15) * 0;
eaB_max(15) = eaB_max(15) * 70;

%% PSwarm
fun = @(x) -my_objective_function3a(x, W, f_ee, r_ee, phi_min, phi_max, t_min, t_max, taskspace_d, radius_tool, radius_scaffold);
Problem = struct('Variables', 15, 'ObjFunction', fun, 'LB', eaB_min, 'UB', eaB_max);

InitPop(1).x = [0; 0; 0; 0; 0; 0; ...
    0; 0; 0; 0; 0; 0; ...
    -length_scaffold; 0; ...
    0];
% InitPop(2).x = [0; 0; 0; 0; 0; 0; ...
%     -dist_tool_b_cg; dist_tool_b_cg ; dist_tool_b_cg; 0; 0; 0; ...
%     -length_scaffold; 0; ...
%     0];
% InitPop(3).x = [0; 120; 240; 0; 120; 240; ...
%     -dist_tool_b_cg; dist_tool_b_cg ; dist_tool_b_cg; 0; 0; 0; ...
%     -length_scaffold; 0; ...
%     0];

% InitPop(1).x = zeros(14,1);
Options = PSwarm('defaults');
Options.Size = 500;
Options.MaxObj = 50000;
Options.MaxIter = 50000;

% Running PSwarm
[BestParticle, BestParticleObj, RunData] = PSwarm(Problem, InitPop, Options);

eaB = BestParticle;
fval = BestParticleObj;

%% Check if Feasible Parameters are found and plotting.
if fval < 0

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

    dist_tooltip = eaB(15);

    draw_cyclops_full(a,B, taskspace, radius_tool, radius_scaffold, length_scaffold, length_overtube, dist_tool_b_cg, dist_tooltip);

    % Workspace Calculation
    f_ee = [0,0,0];
    W = [0,0,0,0,0,0];
    [wp_size, feasible, unfeasible, t] = dex_workspace(a/1000, B/1000, W, f_ee, r_ee/1000, phi_min, phi_max, t_min, t_max);

    plot3(feasible(1,:)*1000, feasible(2,:)*1000, feasible(3,:)*1000, 'b.');
    
    % plot
    figure;
    hold;

    plot3(feasible(1,:), feasible(2,:), feasible(3,:), 'bo');
    plot3(unfeasible(1,:), unfeasible(2,:), unfeasible(3,:), 'r.');

    axis([(min(unfeasible(1,:)) - 20/1000) (max(unfeasible(1,:)) + 20/1000) -(radius_scaffold+5)/1000 (radius_scaffold+5)/1000 -(radius_scaffold+5)/1000 (radius_scaffold+5)/1000]);
else
    X = sprintf('Taskspace not feasible given scaffold and overtube size');
    disp(X);
    
end