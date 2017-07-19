% Given a particular pose, configuration and wrench of the CYCLOPS, determine if
% that pose is feasible
% Pose has 5 components, (x, y, z, alpha, beta)
% Configuration has 2 components, B(2x4), for the base frame attachment
% points and a(2x4), for the position of the attachment points along the
% overtube relative to the centre of mass of the overtube.

function [t, feasible] = feasible_pose(P, a, B, W, f_ee, r_ee, t_min, t_max)

p = P(1:3)';
phi_y = P(4);
phi_z = P(5);


% Rotation Matrix for attachment points
R_y = [cos(phi_y), 0, sin(phi_y); 0, 1, 0; -sin(phi_y), 0, cos(phi_y)];
R_z = [cos(phi_z), -sin(phi_z), 0; sin(phi_z), cos(phi_z), 0; 0, 0, 1];
T_r = R_y * R_z;

% Create Structural Matrix, A
for i=1:size(a, 2)
    
    % Get Attachment point in base frame
    a_b = T_r * a(:,i) + p;
    
    % Find the vector of each tendon
    l = B(:,i) - a_b;
    
    % Get the normalised unit vector
    L(:,i) = l/norm(l);
    
    A(:,i) = [L(:,i); cross(a(:,i),L(:,i))];
end

% Removing DoF for torque in x-direction
% A(4:5,:) = A(5:6,:);
% A(6,:) = [];

% Compute overall wrench
cross_pdt_fee = cross(f_ee', r_ee');
f=[f_ee';cross_pdt_fee]+W';

% % Removing the torque in x-direction
% f(4:5) = f(5:6);
% f(6) = [];

%% Obtaining Tension Solution
% Using analytical method with L1-norm solution
Partition_A = A(:,1:5);
Partition_B = A(:,6);
M = -Partition_A\(f);
N = -Partition_A\Partition_B;

t_low = zeros(6,1);
t_high = zeros(6,1);

for i=1:5
    if N(i) > 0
        t_low(i) = (t_min(i) - M(i))/N(i);
        t_high(i) = (t_max(i) - M(i))/N(i);
    else
        t_low(i) = (t_max(i) - M(i))/N(i);
        t_high(i) = (t_min(i) - M(i))/N(i);
    end
end

t_low(6) = t_min(6);
t_high(6) = t_max(6);

t_B_min = max(t_low);
t_B_max = min(t_high);

if (t_B_min <= t_B_max)
    feasible = 1;
    if t_B_min > t_min(6) && t_B_min < t_max(6)
        t_A = M + N * t_B_min;
        t = [t_A;t_B_min];
    else
        t_A = M + N * t_B_max;
        t = [t_A;t_B_max];
    end
else
    feasible = 0;
    t = zeros(6,1);
end


end

