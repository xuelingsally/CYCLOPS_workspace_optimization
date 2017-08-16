% Spatial rotation matrix around the z axis

function R_z=Rot_z(ang_z)

R_z=[cos(ang_z),-sin(ang_z),0;sin(ang_z),cos(ang_z),0;0,0,1];