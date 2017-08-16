% Calculate roll, pitch, yaw angles from rotation matrix

function [rpy]=rot2rpy(Rot)

% Determining yaw angle
alpha=atan2(Rot(2,1),Rot(1,1));

% Determining pitch angle
beta=atan2(-Rot(3,1),sqrt(Rot(3,2)^2+Rot(3,3)^2));

% Determining roll angle
gamma=atan2(Rot(3,2),Rot(3,3));

% Roll, pitch, yaw vector
rpy=[alpha,beta,gamma];