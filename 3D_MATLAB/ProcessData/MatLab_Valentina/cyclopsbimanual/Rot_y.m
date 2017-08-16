% Spatial rotation matrix around the y axis

function R_y=Rot_y(ang_y)

R_y=[cos(ang_y),0,sin(ang_y);0,1,0;-sin(ang_y),0,cos(ang_y)];