% Internal force vector for Cyclops force control [N]

function [t,W,rank_W]=tension_vector(B,R_des,p_des,r1,r2,r3,r4,r5,r6,r_ee)

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

% Scaling vector
k=10*ones(6,1);

% 6x6 identity matrix
I=eye(6);

% Tendon tension vector
t=pinv(W)*(-f)+(I-pinv(W)*W)*k;


