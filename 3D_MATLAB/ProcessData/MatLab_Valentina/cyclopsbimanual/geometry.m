% Geometry of the cyclops

function [r_tool,s,h1,h2,p_in,T_in,r_ee,p_ee,r1,r2,r3,r4,r5,r5_left,r6,r6_left,B,B_left,P,l,p_x_b]=geometry(L_tool,d_tool,x_f,L,d_bal,d_feed)

% Radius of tool [mm]
r_tool=d_tool/2;

% Distance  between back end of tool and second set of tendon attachment
% points [mm]
x_b=L_tool-L-x_f;

% Radius of the baloon at the tendon feeding point [mm]
r_bal=d_bal/2;

% Distance between P1/P2 and P5 (same for P3/P4 and P6) [mm]
s=sqrt(2)*r_tool;

% Distance between P1 (P2) and P4 (P3) [mm]
h1=sqrt((L^2)+(d_tool^2));

% Distance between P1/P2 and P6 (same for P3/P4 and P5) [mm]
h2=sqrt((L^2)+(2*r_tool^2)); %h2=sqrt((L^2)+(r_tool^2));

% INITIAL CONFIGURATION OF THE TOOL

% Position of the centre of mass of the tool P [mm]
p_in=[0,0,0];

% Orientation of the tool [rad]
R_in=eye(3);

% Transformation matrix between inertial base frame and moving frame at the
% centre of mass of the tool
T_in=[R_in,p_in';zeros(1,3),1];

% Distance between feeding points on baloon and attachment points on tool
% along the x axis of the base frame [mm]
d_x=(d_feed-L)/2; 

% Distance between the centre of mass of the tool and the first set of
% tendon attachment points on the tool along the x axis of the base frame
% [mm]
p_x_f=(L_tool/2)-x_f;

% Distance between the centre of mass of the tool and the second set of
% tendon attachment points on the tool along the x axis of the base frame
% [mm]
p_x_b=(L_tool/2)-x_b;

% Distance between the centre of mass of the tool and the feeding points
% on the baloon (B1,B2,B3,B4) along the z axis of the base frame [mm]
d_z=r_bal;

% Distance between the centre of mass of the tool and the feeding points
% on the baloon (B5, B6) along the y axis of the base frame [mm]
d_y=r_bal;

% Distance between the centre of mass of the tool and the tip of the tool
% [mm]
r_ee=[L_tool/2;0;0];

% Position of the tool end-effector [mm]
p_ee=p_in+(R_in*r_ee)';

% Position of the feeding point on the baloon B1 with respect to the
% inertial base frame [mm]
b1=p_in+[-p_x_f,0,d_z];

% Position of the feeding point on the baloon B2 with respect to the
% inertial base frame [mm]
b2=p_in+[-p_x_f,0,-d_z];

% Position of the feeding point on the baloon B3 with respect to the
% inertial base frame [mm]
b3=p_in+[p_x_b+2*d_x,0,d_z];

% Position of the feeding point on the baloon B4 with respect to the
% inertial base frame [mm]
b4=p_in+[p_x_b+2*d_x,0,-d_z];

% Position of the feeding point on the baloon B5 with respect to the
% inertial base frame [mm]
b5=p_in+[-p_x_f,-d_y,0];

% Position of the feeding point on the baloon B5 left with respect to the
% inertial base frame [mm]
b5_left=p_in+[-p_x_f,d_y,0];

% Position of the feeding point on the baloon B6 with respect to the
% inertial base frame [mm]
b6=p_in+[p_x_b+2*d_x,-d_y,0];

% Position of the feeding point on the baloon B6 left with respect to the
% inertial base frame [mm]
b6_left=p_in+[p_x_b+2*d_x,d_y,0];

% Matrix of feeding point positions on the baloon for the right tool
B=[b1',b2',b3',b4',b5',b6'];

% Matrix of feeding point positions on the baloon for the left tool
B_left=[b1',b2',b3',b4',b5_left',b6_left'];

% Position of the attachment point P1 of the tendon on the tool [mm]
r1=[-p_x_f,0,r_tool]';
p1=p_in'+R_in(1:3,1:3)*r1;

% Position of the attachment point P2 of the tendon on the tool [mm]
r2=[-p_x_f,0,-r_tool]';
p2=p_in'+R_in(1:3,1:3)*r2;

% Position of the attachment point P3 of the tendon on the tool [mm]
r3=[p_x_b,0,r_tool]';
p3=p_in'+R_in(1:3,1:3)*r3;

% Position of the attachment point P4 of the tendon on the tool [mm]
r4=[p_x_b,0,-r_tool]';
p4=p_in'+R_in(1:3,1:3)*r4;

% Position of the attachment point P5 of the tendon on the tool [mm]
r5=[-p_x_f,-r_tool,0]';
p5=p_in'+R_in(1:3,1:3)*r5;

% Position of the attachment point P5 of the tendon on the left tool [mm]
r5_left=[-p_x_f,r_tool,0]';

% Position of the attachment point P6 of the tendon on the tool [mm]
r6=[p_x_b,-r_tool,0]';
p6=p_in'+R_in(1:3,1:3)*r6;

% Position of the attachment point P6 of the tendon on the left tool [mm]
r6_left=[p_x_b,r_tool,0]';

% Matrix of the attachment point positions on the tool
P=[p1,p2,p3,p4,p5,p6];

% Lengths of the tendons
for i=1:6
    l(i)=norm(B(:,i)-P(:,i));
end

