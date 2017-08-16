% Forwards kinematics of the Phantom Omni

function [stylus_pose,T_phantom,T_base,T_link1,T_link2]=Phantom_FK(q,d1,L1,L2)

% DH parameters of the Phantom
d=[d1,0,0];
alpha=[pi/2,0,0];
a=[0,L1,L2];

% Mapping of Phantom joint angles to the kinematics
q_fk(1)=-q(1);
q_fk(2)=q(2);
q_fk(3)=q(3)-q(2)-pi/2;

% Mapping of the Phantom gimbal angles to the kinematics
q_fk(4)=-(q(4)-2.6);
q_fk(5)=q(5)+5.6;
q_fk(6)=-(q(6)-2.6);

% Forward kinematics of the Phantom
T_in=eye(4);
% Position of the stylus
for i=1:3
    T_i=[cos(q_fk(i)),-cos(alpha(i))*sin(q_fk(i)),sin(alpha(i))*sin(q_fk(i)),a(i)*cos(q_fk(i));...
        sin(q_fk(i)),cos(alpha(i))*cos(q_fk(i)),-sin(alpha(i))*cos(q_fk(i)),a(i)*sin(q_fk(i));...
        0,sin(alpha(i)),cos(alpha(i)),d(i);...
        0,0,0,1];
    T=T_in*T_i;
    T_in=T;
    if i==1
        T_base=T_in;
    elseif i==2
        T_link1=T_in;
    elseif i==3
        T_link2=T_in;
    end
end
% Orientation of the stylus
Rot_gimbal=[cos(q_fk(5)),-cos(q_fk(6))*sin(q_fk(5)),sin(q_fk(5))*sin(q_fk(6));
    cos(q_fk(4))*sin(q_fk(5)),cos(q_fk(4))*cos(q_fk(5))*cos(q_fk(6))-sin(q_fk(4))*sin(q_fk(6)),-cos(q_fk(6))*sin(q_fk(4))-cos(q_fk(4))*cos(q_fk(5))*sin(q_fk(6));
    sin(q_fk(4))*sin(q_fk(5)),sin(q_fk(6))*cos(q_fk(4))+sin(q_fk(4))*cos(q_fk(5))*cos(q_fk(6)),cos(q_fk(4))*cos(q_fk(6))-cos(q_fk(5))*sin(q_fk(4))*sin(q_fk(6))];

% Final pose of the stylus
T_phantom=T*[Rot_gimbal,zeros(3,1);zeros(1,3),1];
rpy=rot2rpy(T_phantom(1:3,1:3));
stylus_pose=[T(1,4);T(2,4);T(3,4);rpy'];


