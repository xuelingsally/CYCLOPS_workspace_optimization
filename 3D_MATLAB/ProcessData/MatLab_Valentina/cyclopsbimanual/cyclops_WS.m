% Cyclops controllable workspace

function [A,t_high,t_low,exitflag_high,exitflag_low]=cyclops_WS(t_min,t_max,f_ee,G,B,R_in,p_in,r1,r2,r3,r4,r5,r6,r_ee)

% Allocation of matrix R for computational purposes [m]
R=[r1,r2,r3,r4,r5,r6]/1000;

% Desired position of the tool centre of mass [m]
p_in=p_in/1000;

% Creation of the structure matrix A for the given platform posture
for i=1:6
    % Length vector of each tendon
    u=B(:,i)/1000-p_in-(R_in*R(:,i));
    % Normalized length vector of each tendon
    U(:,i)=u/norm(u);
    % Structure matrix
    A(:,i)=[U(:,i);cross(R(:,i),U(:,i))];
end

% Calculation of the external wrench vector 
% Distance between the centre of mass and the end-effector where the load
% is applied
r_ee=r_ee/1000;
% Resulting wrench vector at the centre of mass
f=[f_ee';cross(f_ee,r_ee)']+G';

% This particular configuration does not allow rotation around the tool
% axis, therefore cannot compensate for applied torques around the x axis
A(4,:)=[];
f(4)=[];

% Initial tension solution
x0=0.1*ones(6,1);

% Optimization function to find the highest allowable tension solution
[t_high,fval_high,exitflag_high]=fmincon(@(x) -myfun2(x),x0,[],[],A,-f,t_min,t_max);

% Optimization function to find the lowest allowable tension solution
[t_low,fval_low,exitflag_low]=fmincon(@(x) myfun2(x),x0,[],[],A,-f,t_min,t_max);
