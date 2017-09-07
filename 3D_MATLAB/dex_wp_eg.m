% Example Usage for dexterous workspace.

a = [3.21	3.21	3.21	-3.21	-3.21	-3.21;
    0.88	-1.75	0.88	0.88	-1.75	0.88;
    1.52	0.00	-1.52	1.52	0.00	-1.52];
  
B = [18.82	18.82	18.82	-18.82	-18.82	-18.82;
     17.50	-35.00	17.50	17.50	-35.00	17.50;
     30.31	0.00	-30.31	30.31	0.00	-30.31];


length_scaffold = max(B(1,:)) - min(B(1,:)); 

f_ee = [0,0,0]; % considered the zero-wrench dexterous workspace
W = [0,0,0,0,0,0];
r_ee = [39.85,0,0];

% Orientation limits
phi_min = [-10/180 * pi; -10/180 * pi];
phi_max = [10/180 * pi; 10/180 * pi];

%Tension limits
t_min = 1 * ones(6,1);
t_max = 60 * ones(6,1);

[dex_wp, feasible, unfeasible,~] = dex_workspace(a/1000, B/1000, W, f_ee, r_ee/1000, phi_min, phi_max, t_min, t_max, length_scaffold/1000);
