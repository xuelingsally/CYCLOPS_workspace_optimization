% Internal force vector for Cyclops force control [N]

function [tau,W,det_N]=internal_force(B,R_in,p_in,r1,r2,r3,r4,r5,r6)

% Allocation of matrix R for computational purposes [m]
R=[r1,r2,r3,r4,r5,r6]/1000;

% Creation of the structure matrix W for the given platform posture
for i=1:6
    % Length vector of each tendon
    u=B(:,i)/1000-p_in-(R_in*R(:,i));
    % Normalized length vector of each tendon
    U(:,i)=u/norm(u);
    % Structure matrix
    W(:,i)=[U(:,i);cross(R(:,i),U(:,i))];
end
% This particular configuration does not allow rotation around the tool
% axis, therefore cannot compensate for applied torques around the x axis
W(4,:)=[];

% The matrix N is a 5x5 square sub-matrix of W
N=W(:,2:6);
det_N=det(N);

% The vector w_last is the last column of W
w_last=W(:,1);

% Set tau_last as a positive value
tau_last=1;

% Internal force vector
tau=[1;-inv(N)*w_last]*tau_last;


