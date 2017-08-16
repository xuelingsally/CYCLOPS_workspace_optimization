% Spatial rotation matrix around the x axis

function R_x=Rot_x(ang_x)

R_x=[1,0,0;0,cos(ang_x),-sin(ang_x);0,sin(ang_x),cos(ang_x)];


