% Inverse kinematics of the cyclops

function [l_diff]=cyclops_IK(alpha,beta,gamma,p_ee,p_ee_new,dr,p_in,T_in,r_ee,r1,r2,r3,r4,r5,r6,B,l)

% Desired rotation matrix
R=rpy2rot(alpha,beta,gamma);

% Desired new position of the tool end-effector [mm]
p_ee_d=p_ee+p_ee_new;

% Desired new position of the tool centre of mass [mm]
p=p_ee_d-(R*r_ee)';

% Transformation matrix between inertial base frame and moving frame at the
% centre of mass of the tool
T=[R,p';zeros(1,3),1];

% Trajectory interpolation
r=[0:dr:1];
r_pt=length(r);
p_des(:,1)=p_in;
p_ee_des(:,1)=p_ee;
for i=1:r_pt-1
    p_des(:,i+1)=p_in*(1-dr*i)+dr*i*p;
    p_ee_des(:,i+1)=p_ee*(1-dr*i)+dr*i*p_ee_d;
end

% Orientation interpolation
Q1=rot2quat(T_in);
Q2=rot2quat(T);
theta = acos(Q1*Q2');
count = 1;
for R_interp=r(:)'
    if theta == 0
        qq = Q1;
    else
        qq = (sin((1-R_interp)*theta) * Q1 + sin(R_interp*theta) * Q2) / sin(theta);
    end
    % Rotation matrix
    R(:,:,count) = quat2rot(qq);
    % Final trajectory
    T_interp(:,:,count) = [R(:,:,count) p_des(:,count); 0 0 0 1];    
    count = count + 1;
end

% Inverse kinematics: calculates the lengths of the tendons for each
% trajectory point followed by the tool tip
for i=1:r_pt
    % Corresponding new position of the tendon attachment points on the tool
    % Position of the attachment point P1 of the tendon on the tool [mm]
    p1_new(:,i)=p_des(:,i)+R(1:3,1:3,i)*r1;
    % Position of the attachment point P2 of the tendon on the tool [mm]
    p2_new(:,i)=p_des(:,i)+R(1:3,1:3,i)*r2;
    % Position of the attachment point P3 of the tendon on the tool [mm]
    p3_new(:,i)=p_des(:,i)+R(1:3,1:3,i)*r3;
    % Position of the attachment point P4 of the tendon on the tool [mm]
    p4_new(:,i)=p_des(:,i)+R(1:3,1:3,i)*r4;
    % Position of the attachment point P5 of the tendon on the tool [mm]
    p5_new(:,i)=p_des(:,i)+R(1:3,1:3,i)*r5;
    % Position of the attachment point P6 of the tendon on the tool [mm]
    p6_new(:,i)=p_des(:,i)+R(1:3,1:3,i)*r6;
    % Matrix of the attachment point positions on the tool
    P_new(:,:,i)=[p1_new(:,i),p2_new(:,i),p3_new(:,i),p4_new(:,i),p5_new(:,i),p6_new(:,i)];
    % Lengths of the tendons
    for j=1:6
        l_new(j,i)=norm(B(:,j)-P_new(:,j,i));
    end
end

% Length difference of the tendons with respect to the homing position
l_diff=l_new-repmat(l',1,r_pt);
