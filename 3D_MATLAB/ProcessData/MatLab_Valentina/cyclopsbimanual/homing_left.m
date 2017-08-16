% Homing configuration of the left tool

function [p_in_left,T_in_left,p_ee_left,P_left,l_left]=homing_left(r_tool,r_ee,r1,r2,r3,r4,r5,r6,B_left)

% Position of the centre of mass of the tool P [mm]
p_in_left=[0,r_tool,0];

% Orientation of the tool [rad]
R_in_left=eye(3);

% Transformation matrix between inertial base frame and moving frame at the
% centre of mass of the tool
T_in_left=[R_in_left,p_in_left';zeros(1,3),1];

% Position of the tool end-effector [mm]
p_ee_left=p_in_left+(R_in_left*r_ee)';

% Position of the attachment point P1 of the tendon on the tool [mm]
p1=p_in_left'+R_in_left(1:3,1:3)*r1;

% Position of the attachment point P2 of the tendon on the tool [mm]
p2=p_in_left'+R_in_left(1:3,1:3)*r2;

% Position of the attachment point P3 of the tendon on the tool [mm]
p3=p_in_left'+R_in_left(1:3,1:3)*r3;

% Position of the attachment point P4 of the tendon on the tool [mm]
p4=p_in_left'+R_in_left(1:3,1:3)*r4;

% Position of the attachment point P5 of the tendon on the tool [mm]
p5=p_in_left'+R_in_left(1:3,1:3)*r5;

% Position of the attachment point P1 of the tendon on the tool [mm]
p6=p_in_left'+R_in_left(1:3,1:3)*r6;

% Matrix of the attachment point positions on the tool
P_left=[p1,p2,p3,p4,p5,p6];

% Lengths of the tendons
for i=1:6
    l_left(i)=norm(B_left(:,i)-P_left(:,i));
end





