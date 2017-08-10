% Internal force vector for Cyclops force control [N]

function [t,W,rank_W]=tension_vector_l1(B,R_des,p_des,r1,r2,r3,r4,r5,r6,r_ee,t_min,t_max)

% Set the force applied at the tool end-effector [N]
f_ee=[0,0,0];

% Set the mass of the tool [kg]
m=0.002;

% Set the acceleration of gravity value [m/(s^2)]
g=9.81;

% Set the gravity force vector [N]
G=[0,0,-m*g,0,0,0];

% Allocation of matrix R for computational purposes [m]
R=[r1,r2,r3,r4,r5,r6]/1000;

% Desired position of the tool centre of mass [m]
p_des=p_des'/1000;

% Calculation of the external wrench vector 
% Distance between the centre of mass and the end-effector where the load
% is applied
r_ee=r_ee/1000;
% Resulting wrench vector at the centre of mass
f=[f_ee';cross(f_ee,r_ee)']+G';

% Creation of the structure matrix W for the given platform posture
for i=1:6
    % Length vector of each tendon
    u=B(:,i)/1000-p_des-(R_des*R(:,i));
    % Normalized length vector of each tendon
    U(:,i)=u/norm(u);
    % Structure matrix
    W(:,i)=[U(:,i);cross(R(:,i),U(:,i))];
end
% This particular configuration does not allow rotation around the tool
% axis, therefore cannot compensate for applied torques around the x axis
W(4,:)=[];
f(4)=[];

% Calculate rank of W to check if the desired configuration is
% singular
rank_W=rank(W);

%% Obtaining Tension Solution
% Using analytical method for L1-norm solution
Partition_A = W(:,1:5);
Partition_B = W(:,6);

% For unfeasible points the matrix approaches singularity
% This line will stop the warning message from MATLAB.
warning('off', 'all');
M = -Partition_A\(f);
N = -Partition_A\Partition_B;
warning('on', 'all');
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
    %feasible = 1;
    if t_B_min > t_min(6) && t_B_min < t_max(6)
        t_A = M + N * t_B_min;
        t = [t_A;t_B_min];
    else
        t_A = M + N * t_B_max;
        t = [t_A;t_B_max];
    end
else
    %feasible = 0;
    t = zeros(6,1);
end


